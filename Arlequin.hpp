//------------------------------------------------------------------------------
// 
//                   Jeferson W D Fernandes and Rodolfo A K Sanches
//                             University of Sao Paulo
//                           (C) 2017 All Rights Reserved
//
// <LicenseText>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//----------------------------------ARLEQUIN------------------------------------
//------------------------------------------------------------------------------

#ifndef ARLEQUIN_H
#define ARLEQUIN_H

#include "Fluid.hpp"
#include "Glue.hpp"

template<int DIM>
class Arlequin{
public:
    typedef Fluid<DIM>                     FluidMesh;
    typedef typename FluidMesh::Elements   Elements;
    typedef typename FluidMesh::Node       Nodes;
    typedef typename FluidMesh::Boundaries Boundary;
    typedef typename Elements::SpecialQuad Quadrature;
    typedef Glue<DIM>                      GlueZone;

private:
    FluidMesh coarseModel, fineModel;

    std::vector<Nodes *>     nodesCoarse_;
    std::vector<Nodes *>     nodesFine_;
    std::vector<Nodes *>     nodesLagrangeFine_;
    std::vector<Nodes *>     nodesLagrangeCoarse_;

    std::vector<Elements *>  elementsCoarse_;
    std::vector<Elements *>  elementsFine_;
    std::vector<GlueZone *>  glueZoneFine_;
    std::vector<GlueZone *>  glueZoneCoarse_;

    std::vector<Boundary *>  boundaryCoarse_;
    std::vector<Boundary *>  boundaryFine_;

    std::vector<int>         elementsGlueZoneFine_;
    std::vector<int>         nodesGlueZoneFine_;
    std::vector<int>         elementsGlueZoneCoarse_;
    std::vector<int>         nodesGlueZoneCoarse_;
    std::vector<int>         elementsFreeZone_;
    std::vector<int>         nodesFreeZone_;

    int numElemCoarse;
    int numElemFine;
    int numElemGlueZoneFine;
    int numElemGlueZoneCoarse;
    int numNodesCoarse;
    int numNodesFine;
    int numNodesGlueZoneFine;
    int numNodesGlueZoneCoarse;
    int numElemFreeZone;
    int numNodesFreeZone;
    int numTimeSteps;
    double dTime;
    int rank;
    int iTimeStep;

    std::pair<idx_t*,idx_t*> domDecompCoarse;//Coarse Model Domain Decomposition
    std::pair<idx_t*,idx_t*> domDecompFine;  //Fine Model Domain Decomposition

    Quadrature quad;

    double pi = M_PI;

public:

    void setFluidModels(FluidMesh coarse, FluidMesh fine);

    int solveSteadyArlequinMovingLaplaceProblem(int iterNumber, 
                                                double tolerance, int steps);

    int solveArlequinProblem(int iterNumber, double tolerance,
                             int problem_type, int time_dependency);
    int solveArlequinProblemMoving(int iterNumber, double tolerance,
                                   int problem_type, int time_dependency);
    int solveArlequinProblemCoarse(int iterNumber, double tolerance,
                                   int problem_type, int time_dependency);

    void setElementBoxes();
    void setCouplingZone();
    void setSignaledDistance();
    void setFreeZone();
    void setWeightFunction(double val);
    double weightFunctionCoarseValue(double r, double epsilon);
    double weightFunctionFineValue(double r, double epsilon);
    void setNodalCorrespondenceFine();
    void setNodalCorrespondenceCoarse();

    std::pair<int,ublas::bounded_vector<double,2> > 
                      searchNodeCorrespondence(typename Nodes::VecLocD x,
                                               std::vector<Nodes *> nodes,
                                               std::vector<Elements *> elements,
                                               int numElem);

    void printVelocity(int step);

};

//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//-------------------------COMPUTE ELEMENT REGIONS/BOXES------------------------
//------------------------------------------------------------------------------
template<>
void Arlequin<2>::setElementBoxes() {
    
    typename Elements::Connectivity connec;
    typename Nodes::VecLocD x1, x2, x3;
    ublas::bounded_vector<double,2> d1, d2, d3, xk, Xk;
    double dCk[3], dck[3];
    std::vector<ublas::bounded_vector<double,2> > di1, di2;

    di1.reserve(3);

    //Compute element boxes for coarse model
    //Only function for straight elements
    for (int jel = 0; jel < numElemCoarse; jel++){
        connec = elementsCoarse_[jel] -> getConnectivity();
        x1 = nodesCoarse_[connec(0)] -> getCoordinates();
        x2 = nodesCoarse_[connec(1)] -> getCoordinates();
        x3 = nodesCoarse_[connec(2)] -> getCoordinates();      

        d1(0) = x2(1) - x1(1);
        d1(1) = x1(0) - x2(0);
        
        di1.push_back(d1);

        d2(0) = x3(1) - x2(1);
        d2(1) = x2(0) - x3(0);

        di1.push_back(d2);

        d3(0) = x1(1) - x3(1);
        d3(1) = x3(0) - x1(0);

        di1.push_back(d3);

        xk(0) = std::min(x1(0),std::min(x2(0), x3(0)));
        xk(1) = std::min(x1(1),std::min(x2(1), x3(1)));

        Xk(0) = std::max(x1(0),std::max(x2(0), x3(0)));
        Xk(1) = std::max(x1(1),std::max(x2(1), x3(1)));        

        dCk[0] = std::max(inner_prod(d1,x1),
                          std::max(inner_prod(d1,x2),inner_prod(d1,x3)));
        dck[0] = std::min(inner_prod(d1,x1),
                          std::min(inner_prod(d1,x2),inner_prod(d1,x3)));

        dCk[1] = std::max(inner_prod(d2,x1),
                          std::max(inner_prod(d2,x2),inner_prod(d2,x3)));
        dck[1] = std::min(inner_prod(d2,x1),
                          std::min(inner_prod(d2,x2),inner_prod(d2,x3)));

        dCk[2] = std::max(inner_prod(d3,x1),
                          std::max(inner_prod(d3,x2),inner_prod(d3,x3)));
        dck[2] = std::min(inner_prod(d3,x1),
                          std::min(inner_prod(d3,x2),inner_prod(d3,x3)));
        
        elementsCoarse_[jel] -> setIntersectionParameters(xk, Xk, dCk, dck,di1);
    };

    di2.reserve(3);
    //Compute element boxes for fine model
    //Only function for straight elements
    for (int jel = 0; jel < numElemFine; jel++){
        connec = elementsFine_[jel] -> getConnectivity();
        x1 = nodesFine_[connec(0)] -> getCoordinates();
        x2 = nodesFine_[connec(1)] -> getCoordinates();
        x3 = nodesFine_[connec(2)] -> getCoordinates();      

        d1(0) = x2(1) - x1(1);
        d1(1) = x1(0) - x2(0);

        di2.push_back(d1);

        d2(0) = x3(1) - x2(1);
        d2(1) = x2(0) - x3(0);
        
        di2.push_back(d2);
        
        d3(0) = x1(1) - x3(1);
        d3(1) = x3(0) - x1(0);

        di2.push_back(d3);

        xk(0) = std::min(x1(0),std::min(x2(0), x3(0)));
        xk(1) = std::min(x1(1),std::min(x2(1), x3(1)));

        Xk(0) = std::max(x1(0),std::max(x2(0), x3(0)));
        Xk(1) = std::max(x1(1),std::max(x2(1), x3(1)));        

        dCk[0] = std::max(inner_prod(d1,x1),
                          std::max(inner_prod(d1,x2),inner_prod(d1,x3)));
        dck[0] = std::min(inner_prod(d1,x1),
                          std::min(inner_prod(d1,x2),inner_prod(d1,x3)));

        dCk[1] = std::max(inner_prod(d2,x1),
                          std::max(inner_prod(d2,x2),inner_prod(d2,x3)));
        dck[1] = std::min(inner_prod(d2,x1),
                          std::min(inner_prod(d2,x2),inner_prod(d2,x3)));

        dCk[2] = std::max(inner_prod(d3,x1),
                          std::max(inner_prod(d3,x2),inner_prod(d3,x3)));
        dck[2] = std::min(inner_prod(d3,x1),
                          std::min(inner_prod(d3,x2),inner_prod(d3,x3)));
        
        elementsFine_[jel] -> setIntersectionParameters(xk, Xk, dCk, dck, di2);
    };

    return;
};


//------------------------------------------------------------------------------
//-------------------COMPUTE NODAL CORRESPONDECE WITH ELEMENTS------------------
//------------------------------------------------------------------------------
template<>
std::pair<int,ublas::bounded_vector<double,2> >
Arlequin<2>::searchNodeCorrespondence(typename Nodes::VecLocD x,
                                      std::vector<Nodes *> nodes, 
                                      std::vector<Elements *> elements, 
                                      int numElem){
    
    typename Nodes::VecLocD xint, x_, deltaX, deltaXsi;
    typename Elements::Connectivity connec;
    ublas::bounded_vector<double,2>            xsi, xsiC;
    QuadShapeFunction<2>                       shapeQuad;
    typename QuadShapeFunction<2>::Values      phi_;
    typename Elements::DimMatrix               ainv;
    double xsiCC[3];
    std::pair<typename Elements::DimVector,typename Elements::DimVector> XK;
    std::pair<double*,double*> dCk;
    std::vector<ublas::bounded_vector<double,2> > di;
    int elemC;
     
    elemC = 150000;
    xsiC(0) = 1.e50;
    xsiC(1) = 1.e50;

    for (int jel = 0; jel < numElem; jel++){
        
        connec = elements[jel] -> getConnectivity();

        //get boxes information        
        XK = elements[jel] -> getXIntersectionParameter();

        //dCk = elements[jel] -> getDIntersectionParameter();
        //di = elements[jel] -> getDiIntersectionParameter();

        //Chech if the node is inside the element box
        if ((x(0) < XK.first(0)) || (x(0) > XK.second(0)) ||
            (x(1) < XK.first(1)) || (x(1) > XK.second(1))) continue;
        
        //Compute nodal correspondence
        xsiCC[0] = 1.e10;
        xsiCC[1] = 1.e10;
        xsiCC[2] = 1.e10;
        
        xsi(0) = 1. / 3.;
        xsi(1) = 1. / 3.;

        shapeQuad.evaluate(xsi,phi_);
        
        x_.clear();
    
        for (int i = 0; i < 6; i++){
            xint = nodes[connec(i)] -> getCoordinates();
            x_(0) += xint(0) * phi_(i);
            x_(1) += xint(1) * phi_(i);                    
        };

        double error = 1.e6;
        
        xsi(0) = 1. / 3.;
        xsi(1) = 1. / 3.;
        
        int iterations = 0;

        while ((error > 1.e-8) && (iterations < 4)) {
            
            iterations++;
            
            deltaX = x - x_;
            
            deltaXsi.clear();
            
            ainv = elements[jel] -> getJacobianMatrixValues(xsi);
        
            noalias(deltaXsi) = prod(trans(ainv),deltaX);
            
            xsi += deltaXsi;
            
            x_.clear();
        
            shapeQuad.evaluate(xsi,phi_);
            
            for (int i=0; i<6; i++){
                xint = nodes[connec(i)] -> getCoordinates();
                x_(0) += xint(0) * phi_(i);
                x_(1) += xint(1) * phi_(i);                    
            };                   
            error = norm_2(deltaXsi);
        };
        
        double t1 = -1.e-2;
        double t2 =  1. - t1;
        
        xsiCC[0] = xsi(0);
        xsiCC[1] = xsi(1);       
        xsiCC[2] = 1. - xsiCC[0] - xsiCC[1];

        if ((xsiCC[0] >= t1) && (xsiCC[1] >= t1) && (xsiCC[2] >= t1) &&
            (xsiCC[0] <= t2) && (xsiCC[1] <= t2) && (xsiCC[2] <= t2)){

            xsiC = xsi;
            elemC = jel;
        };           
    };

    if (fabs(xsi(0)) > 2.) std::cout << "PROBEM SEARCHING NODE CORRESPONDENCE " 
                                     << std::endl;

    return std::make_pair(elemC, xsiC);
};

//------------------------------------------------------------------------------
//--------COMPUTE NODAL CORRESPONDECE OF FINE NODES WITH COARSE ELEMENTS--------
//------------------------------------------------------------------------------
template<>
void Arlequin<2>::setNodalCorrespondenceFine() {

    typename Nodes::VecLocD x;
    ublas::bounded_vector<double,2>            xsi;
    std::pair<int, ublas::bounded_vector<double,2> > corresp;

    //FINE MESH
    for (int inode=0; inode < numNodesFine; inode ++){
        nodesFine_[inode] -> setNodalCorrespondence(1.e10,xsi);
    };

    for (int inode = 0; inode < numNodesGlueZoneFine; inode++) {
        
        x = nodesFine_[nodesGlueZoneFine_[inode]] -> getCoordinates();

        corresp = searchNodeCorrespondence(x, nodesCoarse_, elementsCoarse_, 
                                           elementsCoarse_.size());
   
        nodesFine_[nodesGlueZoneFine_[inode]] -> setNodalCorrespondence(
                                                         corresp.first, 
                                                         corresp.second);
            
        // std::cout << "corresp " << corresp.first 
        //           << " " << corresp.second(0) << " " << corresp.second(1)
        //           << std::endl;

    };

    // for (int i=0; i<numNodesFine; i++){
    //     int elem = nodesFine_[i] -> getNodalElemCorrespondence();
    //     xsi = nodesFine_[i] -> getNodalXsiCorrespondence();
    //     std::cout << "node " << i << " elem " << elem 
    //               << " " << xsi(0) << " " << xsi(1) << std::endl;
    // };

    //Compute correspondence of integration points
    int numberIntPoints = elementsFine_[0] -> getNumberOfIntegrationPoints();
    //if (rank == 0) std::cout << "Int Points " << numberIntPoints << std::endl;

    for (int ielem = 0; ielem < numElemGlueZoneFine; ielem++) {
        
        typename Quadrature::NodalValuesQuad x1, x2;
        typename Nodes::VecLocD x;
        typename Elements::Connectivity connec;
        std::pair<int, ublas::bounded_vector<double,2> > corresp;

        connec = elementsFine_[elementsGlueZoneFine_[ielem]] ->
            getConnectivity();
        
        for (int i = 0; i < 6; i++){
            x = nodesFine_[connec(i)] -> getCoordinates();
            
            x1(i) = x(0);
            x2(i) = x(1);
        };

        int numberIntPoints = elementsFine_[elementsGlueZoneFine_[ielem]] -> 
            getNumberOfIntegrationPoints();

        for (int i = 0; i < numberIntPoints; i++){

            // if (elementsFine_[ielem] -> getIntegPointInGlueZone(i) == true){
                
                x(0) = quad.interpolateQuadraticVariable(x1,i);
                x(1) = quad.interpolateQuadraticVariable(x2,i);
                
            
                corresp = searchNodeCorrespondence(x, nodesCoarse_,
                                                   elementsCoarse_,
                                                   elementsCoarse_.size());
                // if(ielem == 0){
                //     if(i==0)std::cout <<" COOr " << corresp.first << " " << corresp.second(0) << " " << corresp.second(1) << std::endl;
                // };
                elementsFine_[elementsGlueZoneFine_[ielem]] -> 
                    setIntegrationPointCorrespondence(corresp.first,
                                                      corresp.second);
                // };
        };

    };
};

//------------------------------------------------------------------------------
//--------COMPUTE NODAL CORRESPONDECE OF COARSE NODES WITH FINE ELEMENTS--------
//------------------------------------------------------------------------------
template<>
void Arlequin<2>::setNodalCorrespondenceCoarse() {

    typename Nodes::VecLocD x;
    ublas::bounded_vector<double,2>            xsi;
    std::pair<int, ublas::bounded_vector<double,2> > corresp;

    //COARSE MESH
    for (int inode=0; inode < numNodesCoarse; inode ++){
        nodesCoarse_[inode] -> setNodalCorrespondence(1.e10,xsi);
    };

    for (int inode = 0; inode < numNodesGlueZoneFine; inode++) {
        
        x = nodesCoarse_[nodesGlueZoneFine_[inode]] -> getCoordinates();
        
        corresp = searchNodeCorrespondence(x, nodesFine_, elementsFine_, 
                                           elementsFine_.size());
   
        nodesCoarse_[nodesGlueZoneFine_[inode]] -> setNodalCorrespondence(
                                                           corresp.first, 
                                                           corresp.second);
    };

    // for (int i=0; i<numNodesFine; i++){
    //     int elem = nodesFine_[i] -> getNodalElemCorrespondence();
    //     xsi = nodesFine_[i] -> getNodalXsiCorrespondence();
    //     std::cout << "node " << i << " elem " << elem 
    //               << " " << xsi(0) << " " << xsi(1) << std::endl;
    // };

    //Compute correspondence of integration points
    int numberIntPoints = elementsCoarse_[0] -> getNumberOfIntegrationPoints();
    std::cout << "inT Points " << numberIntPoints << std::endl;

    for (int ielem = 0; ielem < numElemGlueZoneFine; ielem++) {
        
        typename Quadrature::NodalValuesQuad x1, x2;
        typename Nodes::VecLocD x;
        typename Elements::Connectivity connec;
        std::pair<int, ublas::bounded_vector<double,2> > corresp;

        connec = elementsCoarse_[elementsGlueZoneFine_[ielem]] -> 
            getConnectivity();
        
        for (int i = 0; i < 6; i++){
            x = nodesCoarse_[connec(i)] -> getCoordinates();
            
            x1(i) = x(0);
            x2(i) = x(1);
        };

        int numberIntPoints = elementsCoarse_[elementsGlueZoneFine_[ielem]] -> 
            getNumberOfIntegrationPoints();

        for (int i = 0; i < numberIntPoints; i++){

            if (elementsCoarse_[ielem] -> 
                getIntegPointInGlueZone(i) == true){
                
                x(0) = quad.interpolateQuadraticVariable(x1,i);
                x(1) = quad.interpolateQuadraticVariable(x2,i);
                
                corresp = searchNodeCorrespondence(x, nodesFine_, elementsFine_,
                                                   elementsFine_.size());
                
                
                elementsCoarse_[elementsGlueZoneFine_[ielem]] -> 
                    setIntegrationPointCorrespondence(corresp.first,
                                                      corresp.second);
            };
        };
        // std::cout << "corresp " << corresp.first 
        //<< " " << corresp.second(0) << std::endl;
    };
};

//------------------------------------------------------------------------------
//--------------------SETS THE COUPLING ZONE IN COARSE MODEL--------------------
//------------------------------------------------------------------------------
template<>
void Arlequin<2>::setSignaledDistance(){

    typename Elements::Connectivity connec;
    typename Nodes::VecLocD x;

    // Set Dist Function
    for (int ino = 0; ino < numNodesFine; ino++){
        x = nodesFine_[ino] -> getCoordinates();
        
        double swd = 3.;
        double cswd = 8.;//1.-swd;
        double distx1 = -x(0)+swd;
        double distx2 = -cswd + x(0);
        double disty1 = -x(1)+swd;
        double disty2 = -cswd + x(1);
        
        double dist = 0.;
        if((x(0)<=swd) && (x(1)<=swd)){
            dist = sqrt((x(0)-swd)*(x(0)-swd)+(x(1)-swd)*(x(1)-swd));
        };
        if((x(0)<=swd) && (x(1)>=cswd)){
            dist = sqrt((x(0)-swd)*(x(0)-swd)+(x(1)-cswd)*(x(1)-cswd));
        };
        if((x(0)>=cswd) && (x(1)>=cswd)){
            dist = sqrt((x(0)-cswd)*(x(0)-cswd)+(x(1)-cswd)*(x(1)-cswd));
        };
        if((x(0)>=cswd) && (x(1)<=swd)){
            dist = sqrt((x(0)-cswd)*(x(0)-cswd)+(x(1)-swd)*(x(1)-swd));
        };
        if((x(0)<=swd) && ((x(1)>swd)&&(x(1)<cswd))){
            dist = distx1;
        };
        if((x(0)>=cswd) && ((x(1)>swd)&&(x(1)<cswd))){
            dist = distx2;
        };        
        if((x(1)<=swd) && ((x(0)>swd)&&(x(0)<cswd))){
            dist = disty1;
        };
        if((x(1)>=cswd) && ((x(0)>swd)&&(x(0)<cswd))){
            dist = disty2;
        };        
        if(((x(1)>swd)&&(x(1)<cswd)) && ((x(0)>swd)&&(x(0)<cswd))){
            dist = distx1;
            if (fabs(dist)>=fabs(distx2)){
                dist=distx2;
            };
            if (fabs(dist)>=fabs(disty1)){
                dist=disty1;
            };
            if (fabs(dist)>=fabs(disty2)){
                dist=disty2;
            };
        };
        nodesFine_[ino] -> setDistFunction(dist);
    };    

    for (int i = 0; i < numElemFine; i++){
        for (int j = 0; 
             j < elementsFine_[i] -> getNumberOfIntegrationPoints(); j++){
        
            x = elementsFine_[i] -> getIntegPointCoordinatesValue(j);
        
            double swd = 3.;
            double cswd = 8.;//1.-swd;
            double distx1 = -x(0)+swd;
            double distx2 = -cswd + x(0);
            double disty1 = -x(1)+swd;
            double disty2 = -cswd + x(1);
            
            double dist = 0.;
            if((x(0)<=swd) && (x(1)<=swd)){
                dist = sqrt((x(0)-swd)*(x(0)-swd)+(x(1)-swd)*(x(1)-swd));
            };
            if((x(0)<=swd) && (x(1)>=cswd)){
                dist = sqrt((x(0)-swd)*(x(0)-swd)+(x(1)-cswd)*(x(1)-cswd));
            };
            if((x(0)>=cswd) && (x(1)>=cswd)){
                dist = sqrt((x(0)-cswd)*(x(0)-cswd)+(x(1)-cswd)*(x(1)-cswd));
            };
            if((x(0)>=cswd) && (x(1)<=swd)){
                dist = sqrt((x(0)-cswd)*(x(0)-cswd)+(x(1)-swd)*(x(1)-swd));
            };
            if((x(0)<=swd) && ((x(1)>swd)&&(x(1)<cswd))){
                dist = distx1;
            };
            if((x(0)>=cswd) && ((x(1)>swd)&&(x(1)<cswd))){
                dist = distx2;
            };        
            if((x(1)<=swd) && ((x(0)>swd)&&(x(0)<cswd))){
                dist = disty1;
            };
            if((x(1)>=cswd) && ((x(0)>swd)&&(x(0)<cswd))){
                dist = disty2;
            };        
            if(((x(1)>swd)&&(x(1)<cswd)) && ((x(0)>swd)&&(x(0)<cswd))){
                dist = distx1;
                if (fabs(dist)>=fabs(distx2)){
                    dist=distx2;
                };
                if (fabs(dist)>=fabs(disty1)){
                    dist=disty1;
                };
                if (fabs(dist)>=fabs(disty2)){
                    dist=disty2;
                };
            };
            elementsFine_[i] -> setIntegPointDistFunction(j,dist);
        };    
    };




    // Coarse
    for (int ino = 0; ino < numNodesCoarse; ino++){
        x = nodesCoarse_[ino] -> getCoordinates();

        double h = .5;
        double T = 50.;        
        double swd = 3.;
        double cswd = 8.;//1.-swd;
        double swdy = 3. + h * sin(pi * iTimeStep * dTime / T);
        double cswdy = 8. + h * sin(pi * iTimeStep * dTime / T);//1.-swd;
        double distx1 = -x(0)+swd;
        double distx2 = -cswd + x(0);
        double disty1 = -x(1)+(swdy);
        double disty2 = -(cswdy) + x(1);
        
        double dist = 0.;
        if((x(0)<=swd) && (x(1)<=swdy)){
            dist = sqrt((x(0)-swd)*(x(0)-swd)+(x(1)-swdy)*(x(1)-swdy));
        };
        if((x(0)<=swd) && (x(1)>=cswdy)){
            dist = sqrt((x(0)-swd)*(x(0)-swd)+(x(1)-cswdy)*(x(1)-cswdy));
        };
        if((x(0)>=cswd) && (x(1)>=cswdy)){
            dist = sqrt((x(0)-cswd)*(x(0)-cswd)+(x(1)-cswdy)*(x(1)-cswdy));
        };
        if((x(0)>=cswd) && (x(1)<=swdy)){
            dist = sqrt((x(0)-cswd)*(x(0)-cswd)+(x(1)-swdy)*(x(1)-swdy));
        };
        if((x(0)<=swd) && ((x(1)>swdy)&&(x(1)<cswdy))){
            dist = distx1;
        };
        if((x(0)>=cswd) && ((x(1)>swdy)&&(x(1)<cswdy))){
            dist = distx2;
        };        
        if((x(1)<=swdy) && ((x(0)>swd)&&(x(0)<cswd))){
            dist = disty1;
        };
        if((x(1)>=cswdy) && ((x(0)>swd)&&(x(0)<cswd))){
            dist = disty2;
        };        
        if(((x(1)>swdy)&&(x(1)<cswdy)) && ((x(0)>swd)&&(x(0)<cswd))){
            dist = distx1;
            if (fabs(dist)>=fabs(distx2)){
                dist=distx2;
            };
            if (fabs(dist)>=fabs(disty1)){
                dist=disty1;
            };
            if (fabs(dist)>=fabs(disty2)){
                dist=disty2;
            };
        };
        nodesCoarse_[ino] -> setDistFunction(dist);
    };  

    for (int i = 0; i < numElemCoarse; i++){
        for (int j = 0; 
             j < elementsCoarse_[i] -> getNumberOfIntegrationPoints(); j++){
        
            x = elementsCoarse_[i] -> getIntegPointCoordinatesValue(j);
        
            double h = .5;
            double T = 50.;    

            double swd = 3.;
            double cswd = 8.;//1.-swd;
            double swdy = 3. + h * sin(pi * iTimeStep * dTime / T);
            double cswdy = 8. + h * sin(pi * iTimeStep * dTime / T);//1.-swd;
            double distx1 = -x(0)+swd;
            double distx2 = -cswd + x(0);
            double disty1 = -x(1)+(swd);
            double disty2 = -(cswd) + x(1);
            
            double dist = 0.;
            if((x(0)<=swd) && (x(1)<=swdy)){
                dist = sqrt((x(0)-swd)*(x(0)-swd)+(x(1)-swdy)*(x(1)-swdy));
            };
            if((x(0)<=swd) && (x(1)>=cswdy)){
                dist = sqrt((x(0)-swd)*(x(0)-swd)+(x(1)-cswdy)*(x(1)-cswdy));
            };
            if((x(0)>=cswd) && (x(1)>=cswdy)){
                dist = sqrt((x(0)-cswd)*(x(0)-cswd)+(x(1)-cswdy)*(x(1)-cswdy));
            };
            if((x(0)>=cswd) && (x(1)<=swdy)){
                dist = sqrt((x(0)-cswd)*(x(0)-cswd)+(x(1)-swdy)*(x(1)-swdy));
            };
            if((x(0)<=swd) && ((x(1)>swdy)&&(x(1)<cswdy))){
                dist = distx1;
            };
            if((x(0)>=cswd) && ((x(1)>swdy)&&(x(1)<cswdy))){
                dist = distx2;
            };        
            if((x(1)<=swdy) && ((x(0)>swd)&&(x(0)<cswd))){
                dist = disty1;
            };
            if((x(1)>=cswdy) && ((x(0)>swd)&&(x(0)<cswd))){
                dist = disty2;
            };        
            if(((x(1)>swdy)&&(x(1)<cswdy)) && ((x(0)>swd)&&(x(0)<cswd))){
                dist = distx1;
                if (fabs(dist)>=fabs(distx2)){
                    dist=distx2;
                };
                if (fabs(dist)>=fabs(disty1)){
                    dist=disty1;
                };
                if (fabs(dist)>=fabs(disty2)){
                    dist=disty2;
                };
            };
            elementsCoarse_[i] -> setIntegPointDistFunction(j,dist);
        };    
    }; 

};

//------------------------------------------------------------------------------
//--------------------SETS THE COUPLING ZONE IN COARSE MODEL--------------------
//------------------------------------------------------------------------------
template<>
void Arlequin<2>::setCouplingZone(){

    typename Elements::Connectivity connec;
    typename Nodes::VecLocD x;
    // double dist;
    int flag;
    int nodesCZ[numNodesFine];
    int nodesCZ2[numNodesFine];

    double lim1 = 2.99;
    double lim2 = 8.01;
    double tick = 0.99;

    for (int i = 0; i < numNodesFine; i++) nodesCZ[i] = 0;    

    elementsGlueZoneFine_.reserve(numElemFine / 3);
    nodesGlueZoneFine_.reserve(numNodesFine / 3);
    glueZoneFine_.reserve(numNodesFine / 3);
    int index = 0;
  
    //Defines a criterion to select the elements that are in the glue zone
    for (int jel = 0; jel < numElemFine; jel++){
        
        connec = elementsFine_[jel] -> getConnectivity();
        flag = 0;

        for (int ino = 0; ino < 6; ino++){
            x = nodesFine_[connec(ino)] -> getCoordinates();

            if ((x(0) < lim1) || (x(0) > lim2) || 
                (x(1) < lim1) || (x(1) > lim2)){
                flag = 1;
                break;
            };
        };

        if (flag > 0) {
            elementsGlueZoneFine_.push_back(jel);
            elementsFine_[jel] -> setGlueZone();

            GlueZone *el = new GlueZone(index++,jel);
            glueZoneFine_.push_back(el);
            
            for (int i=0; 
                 i < elementsFine_[jel] -> getNumberOfIntegrationPoints(); i++){
                
                x = elementsFine_[jel] -> getIntegPointCoordinatesValue(i);  

                if ((x(0) < lim1) || (x(0) > lim2) || 
                    (x(1) < lim1) || (x(1) > lim2)){
                    elementsFine_[jel] -> setIntegPointInGlueZone(i);    
                };
            };
        };        
    };

  
    //Defines which nodes are in the glue zone
    numElemGlueZoneFine = elementsGlueZoneFine_.size();
    for (int i = 0; i < numElemGlueZoneFine; i++){
        connec = elementsFine_[elementsGlueZoneFine_[i]] -> getConnectivity();

        for (int ino = 0; ino < 6; ino++){
            nodesCZ[connec(ino)] += 1;
        };
    };

    //Compute number of nodes in the glue zone
    numNodesGlueZoneFine = 0;
    for (int i = 0; i < numNodesFine; i++){
        if(nodesCZ[i] > 0) {
            numNodesGlueZoneFine += 1;
            nodesGlueZoneFine_.push_back(i);
        };
    };
      
    if (rank == 0) std::cout << "GLUE ZONE - Nodes = " << numNodesGlueZoneFine 
                             << " Elements = " 
                             << elementsGlueZoneFine_.size() << std::endl;
    
    for (int i = 0; i < numNodesGlueZoneFine; i++){
        typename Nodes::VecLocD x;
        x = nodesFine_[nodesGlueZoneFine_[i]] -> getCoordinates();
        
        Nodes *no = new Nodes(x,i);
        nodesLagrangeFine_.push_back(no);
    };

    for (int i = 0; i < numElemGlueZoneFine; i++){
        typename Elements::Connectivity connecAux;

        connec = elementsFine_[elementsGlueZoneFine_[i]] -> getConnectivity();
        
        for (int ino = 0; ino < numNodesGlueZoneFine; ino++){
            if (nodesGlueZoneFine_[ino] == connec(0)) connecAux(0) = ino;
            if (nodesGlueZoneFine_[ino] == connec(1)) connecAux(1) = ino;
            if (nodesGlueZoneFine_[ino] == connec(2)) connecAux(2) = ino;
            if (nodesGlueZoneFine_[ino] == connec(3)) connecAux(3) = ino;
            if (nodesGlueZoneFine_[ino] == connec(4)) connecAux(4) = ino;
            if (nodesGlueZoneFine_[ino] == connec(5)) connecAux(5) = ino;
        };
        
        glueZoneFine_[i] -> setConnectivity(connecAux);
        glueZoneFine_[i] -> setNodes(nodesLagrangeFine_);

    };
    std::cout << "AQUI 1 " << rank << std::endl;     
    MPI_Barrier(PETSC_COMM_WORLD);      

    // Glue Zone in coarse mesh

    for (int i = 0; i < numNodesCoarse; i++) nodesCZ2[i] = 0;    

    elementsGlueZoneCoarse_.reserve(numElemCoarse / 3);
    nodesGlueZoneCoarse_.reserve(numNodesCoarse / 3);

    std::cout << "AQUI 2 " << rank << std::endl;     
    MPI_Barrier(PETSC_COMM_WORLD);          

    //Defines a criterion to select the elements that are in the glue zone
    // for (int jel = 0; jel < numElemCoarse; jel++){
        
    //     connec = elementsCoarse_[jel] -> getConnectivity();
    //     flag = 0;

    //     for (int ino = 0; ino < 6; ino++){
    //         x = nodesCoarse_[connec(ino)] -> getCoordinates();

    //         if ((x(0) < lim1) || (x(0) > lim2) || 
    //             (x(1) < lim1) || (x(1) > lim2)){

    //         }else{
    //             if ((x(0) > lim1 + tick) && (x(0) < lim2 - tick) &&
    //                 (x(1) > lim1 + tick) && (x(1) < lim2 - tick)){

    //             }else{
    //                 flag = 1;
    //                 break;
    //             };
    //         };
    //     };
    //     if (flag > 0) {
    //         elementsGlueZoneCoarse_.push_back(jel);
    //         elementsCoarse_[jel] -> setGlueZone();

    //         GlueZone *el = new GlueZone(index++,jel);
    //         glueZoneCoarse_.push_back(el);
            
    //         for (int i=0; i < elementsCoarse_[jel] -> 
    //                  getNumberOfIntegrationPoints(); i++){
                
    //             x = elementsCoarse_[jel] -> getIntegPointCoordinatesValue(i);  

    //             if ((x(0) < lim1) || (x(0) > lim2) || 
    //                 (x(1) < lim1) || (x(1) > lim2)){

    //             }else{
    //                 if ((x(0) > lim1 + tick) && (x(0) < lim2 - tick) &&
    //                     (x(1) > lim1 + tick) && (x(1) < lim2 - tick)){

    //                 }else{
    //                 elementsCoarse_[jel] -> setIntegPointInGlueZone(i);
    //                 };
    //             };
    //         };
    //     };        
    // };
    // std::cout << "AQUI 3 " << rank << std::endl;     
    // MPI_Barrier(PETSC_COMM_WORLD);      
    // //Defines which nodes are in the glue zone
    // numElemGlueZoneCoarse = elementsGlueZoneCoarse_.size();
    // for (int i = 0; i < numElemGlueZoneCoarse; i++){
    //     connec = elementsCoarse_[elementsGlueZoneCoarse_[i]] -> 
    //         getConnectivity();

    //     for (int ino = 0; ino < 6; ino++){
    //         nodesCZ2[connec(ino)] += 1;
    //     };
    // };
    // std::cout << "AQUI 3.5 " << rank << std::endl;     
    // MPI_Barrier(PETSC_COMM_WORLD);      

    // //Compute number of nodes in the glue zone
    // numNodesGlueZoneCoarse = 0;
    // for (int i = 0; i < numNodesCoarse; i++){
    //     if(nodesCZ2[i] > 0) {
    //         numNodesGlueZoneCoarse += 1;
    //         nodesGlueZoneCoarse_.push_back(i);
    //     };
    // };
    //     std::cout << "AQUI 4 " << rank << std::endl;     
    // MPI_Barrier(PETSC_COMM_WORLD);      
    // if (rank == 0) std::cout << "GLUE ZONE - Nodes = " << numNodesGlueZoneCoarse
    //                          << " Elements = " 
    //                          << elementsGlueZoneCoarse_.size() << std::endl;

    // for (int i = 0; i < numNodesGlueZoneCoarse; i++){
    //     typename Nodes::VecLocD x;
    //     x = nodesCoarse_[nodesGlueZoneCoarse_[i]] -> getCoordinates();
        
    //     Nodes *no = new Nodes(x,i);
    //     nodesLagrangeCoarse_.push_back(no);
    // };

    // for (int i = 0; i < numElemGlueZoneCoarse; i++){
    //     typename Elements::Connectivity connecAux;

    //     connec = elementsCoarse_[elementsGlueZoneCoarse_[i]] -> 
    //         getConnectivity();
        
    //     for (int ino = 0; ino < numNodesGlueZoneCoarse; ino++){
    //         if (nodesGlueZoneCoarse_[ino] == connec(0)) connecAux(0) = ino;
    //         if (nodesGlueZoneCoarse_[ino] == connec(1)) connecAux(1) = ino;
    //         if (nodesGlueZoneCoarse_[ino] == connec(2)) connecAux(2) = ino;
    //         if (nodesGlueZoneCoarse_[ino] == connec(3)) connecAux(3) = ino;
    //         if (nodesGlueZoneCoarse_[ino] == connec(4)) connecAux(4) = ino;
    //         if (nodesGlueZoneCoarse_[ino] == connec(5)) connecAux(5) = ino;
    //     };
        
    //     glueZoneCoarse_[i] -> setConnectivity(connecAux);
    //     glueZoneCoarse_[i] -> setNodes(nodesLagrangeCoarse_);

    // };
     

    setSignaledDistance();

    std::cout << "AQUI 6 " << rank << std::endl;     
    MPI_Barrier(PETSC_COMM_WORLD);      
    std::cout << "AQUI 7 " << rank << std::endl;     
};


//------------------------------------------------------------------------------
//----------------------SETS THE FREE ZONE IN COARSE- MODEL---------------------
//------------------------------------------------------------------------------
template<>
void Arlequin<2>::setFreeZone(){

    typename Elements::Connectivity connec;
    typename Nodes::VecLocD x;
    // double dist;
    int flag;

    elementsFreeZone_.clear();
    nodesFreeZone_.clear();

    elementsFreeZone_.reserve(numElemCoarse / 5);
    nodesFreeZone_.reserve(numNodesCoarse / 5);

    double h = .5;
    double T = 50.;
    double d1 = 3.;
    double d2 = 8.;
    double d3 = 3. + h * sin(pi * iTimeStep * dTime / T);
    double d4 = 8. + h * sin(pi * iTimeStep * dTime / T);

    //Defines a criterion to select the elements that are in the glue zone
    for (int jel = 0; jel < numElemCoarse; jel++){
                
        elementsCoarse_[jel] -> clearCompressibility();

        connec = elementsCoarse_[jel] -> getConnectivity();
        flag = 0;

        for (int ino = 0; ino < 6; ino++){
            x = nodesCoarse_[connec(ino)] -> getCoordinates();

            if ((x(0) >= d1) && (x(0) <= d2) &&
                (x(1) >= d3) && (x(1) <= d4))  {
                flag += 1;
            };
        };

        if (flag == 6) {
            elementsFreeZone_.push_back(jel);
            //elementsCoarse_[jel] -> setCompressibility();
            //std::cout << "Element " << jel << std::endl;
        }; 
        if (flag == 6) {
            nodesFreeZone_.push_back(connec(0));
            nodesFreeZone_.push_back(connec(1));
            nodesFreeZone_.push_back(connec(2));
            nodesFreeZone_.push_back(connec(3));
            nodesFreeZone_.push_back(connec(4));
            nodesFreeZone_.push_back(connec(5));
        };

    };

    // for (int jel = 0; jel < numElemCoarse; jel++){
                
    //     connec = elementsCoarse_[jel] -> getConnectivity();
    //     flag = 0;
        
    //     for (int ino = 0; ino < 6; ino++){
    //         x = nodesCoarse_[connec(ino)] -> getCoordinates();
            
    //         if ((x(0) < d1) || (x(0) > d2) ||
    //             (x(1) < d1) || (x(1) > d2))  {
    //             flag++;
    //         };
    //     };
    //     if (flag == 6
    //     nodesFreeZone_.push_back(ino);
    // };

    numElemFreeZone = elementsFreeZone_.size();
    numNodesFreeZone = nodesFreeZone_.size();
};

//------------------------------------------------------------------------------
//----------------------COMPUTES THE WEIGHT FUNCTION VALUE----------------------
//------------------------------------------------------------------------------
template<>
double Arlequin<2>::weightFunctionFineValue(double r, double epsilon){

    double wFuncValue;

    if (r <= 2.5){
        wFuncValue = 1. - epsilon;
    } else {
        if (r > 3.55){
            wFuncValue = 0.;
        } else {
            wFuncValue = 1. - epsilon;

            //Linear weight function
            // 0.5 < r < 3.5
            //wFuncValue = - (1. - epsilon) * (r - 3.5) / 3.;
            // 2.5 < r < 3.5
            //wFuncValue = - (1. - epsilon) * (r - 3.5);
            // 2.0 < r < 3.5
            //wFuncValue = -(1 - epsilon) * (r - 3.5) / 1.5;
            // 1.5 < r < 3.5
            //wFuncValue = -(1. - epsilon) * (r - 3.5) / 2.;
            // 3.0 < r < 3.5
            //wFuncValue = 2. * (1. - epsilon) * (r - 3.5) + 1.;
            // 1.5 < r < 3.0
            //wFuncValue = -(1. - epsilon) * (r - 3) / 1.5;
            // 2.0 < r < 3.0
            //wFuncValue = -(1. - epsilon) * (r - 3) / 2.0;
            //wFuncValue = 1. - epsilon;

            //Cubic weight function
            // 2.5 < r < 3.5
            // wFuncValue = 1. - 2. * (epsilon - 1.) * r * r * r 
            //      - 18. * (1. - epsilon) * r * r
            //      - 52.5 * (epsilon - 1.) * r
            //      - 50. + 49. * epsilon;
            // 1.5 < r < 3.5
            // wFuncValue = 1. - 0.25 * (epsilon - 1.) * r * r * r
            //     - 15. / 8. * (1. - epsilon) * r * r
            //     - 63. / 16. * (epsilon - 1.) * r
            //     - (81. - 49. * epsilon) / 32.;
            // 2.0 < r < 3.5
            // wFuncValue = 1. - 16. / 27. * (epsilon - 1.) * r * r * r
            //     - 44. / 9. * (1. - epsilon) * r * r 
            //     - 112. / 9. * (epsilon - 1.) * r
            //     - (272. - 245. * epsilon) / 27.;
            // 3.0 < r < 3.5
            // wFuncValue = 1. - 16. * (epsilon - 1.) * r * r * r 
            //     - 156. * (1. - epsilon) * r * r 
            //     - 504. * (epsilon - 1.) * r
            //     - (540. - 539. * epsilon);


            //Quadratic weight function
            // 2.5 < r < 3.5
            // if (r <= 3.){
            //     wFuncValue = 1. - 2. * (1. - epsilon) * r * r
            //         - 10. * (epsilon - 1.) * r
            //         - 0.5 * (25. - 23. * epsilon);
            // }else{
            //     wFuncValue = 1. - 2. * (epsilon - 1.) * r * r
            //         - 14. * (1. - epsilon) * r
            //         - 0.5 * (49. * epsilon - 47.);
            // };
            // 3.0 < r < 3.5
            // if (r <= 3.25){
            //     wFuncValue = 1. - 8. * (1. - epsilon) * r * r
            //         - 48. * (epsilon - 1.) * r
            //         - (72. - 71. * epsilon);
            // }else{
            //     wFuncValue = 1. - 8. * (epsilon - 1.) * r * r
            //         - 56. * (1. - epsilon) * r
            //         - (98. * epsilon - 97.);
            // };
            // // 2.0 < r < 3.5
            // if (r <= 2.75){
            //     wFuncValue = 1. - 8. / 9. * (1. - epsilon) * r * r
            //         - 32. / 9. * (epsilon - 1.) * r
            //         - (32. - 23. * epsilon) / 9.;
            // }else{
            //     wFuncValue = 1. - 8. / 9. * (epsilon - 1.) * r * r
            //         - 56. / 9. * (1. - epsilon) * r
            //         - (98. * epsilon - 89.) / 9.;
            // };
            // 1.5 < r < 3.5
            // if (r <= 2.5){
            //     wFuncValue = 1. - 1. / 2. * (1. - epsilon) * r * r
            //         - 3. / 2. * (epsilon - 1.) * r
            //         - (1.125 - 0.125 * epsilon);
            // }else{
            //     wFuncValue = 1. - 1. / 2. * (epsilon - 1.) * r * r
            //         - 7. / 2. * (1. - epsilon) * r
            //         - (6.125 * epsilon - 5.125);
            // };
        };
    };

    return wFuncValue;
        
};

template<>
double Arlequin<2>::weightFunctionCoarseValue(double r, double epsilon){

    double wFuncValue;

    if (r <= 2.5){
        wFuncValue = epsilon;
    } else {
        if (r >= 3.5){
            wFuncValue = 1.;
        } else {
            wFuncValue = epsilon;

            //Linear weight function
            // 0.5 < r < 3.5
            //wFuncValue = (1. - epsilon) * (r - 3.5) / 3. + 1.;
            // 2.5 < r < 3.5
            //wFuncValue = (1. - epsilon) * (r - 3.5) + 1.;
            // 2.0 < r < 3.5
            //wFuncValue = (1 - epsilon) * (r - 3.5) / 1.5 + 1.;
            // 1.5 < r < 3.5
            //wFuncValue = (1. - epsilon) * (r - 3.5) / 2. + 1.;
            // 3.0 < r < 3.5
            //wFuncValue = 2. * (1. - epsilon) * (r - 3.5) + 1.;
            // 1.5 < r < 3.0
            //wFuncValue = (1. - epsilon) * (r - 3) / 1.5 + 1.;
            // 2.0 < r < 3.0
            //wFuncValue = (1. - epsilon) * (r - 3) / 2.0 + 1.;
            //wFuncValue = epsilon;

            //Cubic weight function
            // 2.5 < r < 3.5
            // wFuncValue = 2. * (epsilon - 1.) * r * r * r 
            //     + 18. * (1. - epsilon) * r * r
            //     + 52.5 * (epsilon - 1.) * r
            //     + 50. - 49. * epsilon;
            // 1.5 < r < 3.5
            // wFuncValue = 0.25 * (epsilon - 1.) * r * r * r
            //     + 15. / 8. * (1. - epsilon) * r * r
            //     + 63. / 16. * (epsilon - 1.) * r
            //     + (81. - 49. * epsilon) / 32.;see
            // 2.0 < r < 3.5
            // wFuncValue = 16. / 27. * (epsilon - 1.) * r * r * r
            //     + 44. / 9. * (1. - epsilon) * r * r 
            //     + 112. / 9. * (epsilon - 1.) * r
            //     + (272. - 245. * epsilon) / 27.;
            // 3.0 < r < 3.5
            // wFuncValue = 16. * (epsilon - 1.) * r * r * r 
            //     + 156. * (1. - epsilon) * r * r 
            //     + 504. * (epsilon - 1.) * r
            //     + (540. - 539. * epsilon);

            //Quadratic weight function
            // 2.5 < r < 3.5
            // if (r <= 3.){
            //     wFuncValue = 2. * (1. - epsilon) * r * r
            //         + 10. * (epsilon - 1.) * r
            //         + 0.5 * (25. - 23. * epsilon);
            // }else{
            //     wFuncValue = 2. * (epsilon - 1.) * r * r
            //         + 14. * (1. - epsilon) * r
            //         + 0.5 * (49. * epsilon - 47.);
            // };
            // 3.0 < r < 3.5
            // if (r <= 3.25){
            //     wFuncValue = 8. * (1. - epsilon) * r * r
            //         + 48. * (epsilon - 1.) * r
            //         + (72. - 71. * epsilon);
            // }else{
            //     wFuncValue = 8. * (epsilon - 1.) * r * r
            //         + 56. * (1. - epsilon) * r
            //         + (98. * epsilon - 97.);
            // };
            // // 2.0 < r < 3.5
            // if (r <= 2.75){
            //     wFuncValue = 8. / 9. * (1. - epsilon) * r * r
            //         + 32. / 9. * (epsilon - 1.) * r
            //         + (32. - 23. * epsilon) / 9.;
            // }else{
            //     wFuncValue = 8. / 9. * (epsilon - 1.) * r * r
            //         + 56. / 9. * (1. - epsilon) * r
            //         + (98. * epsilon - 89.) / 9.;
            // };
            // 1.5 < r < 3.5
            // if (r <= 2.5){
            //     wFuncValue = 1. / 2. * (1. - epsilon) * r * r
            //         + 3. / 2. * (epsilon - 1.) * r
            //         + (1.125 - 0.125 * epsilon);
            // }else{
            //     wFuncValue = 1. / 2. * (epsilon - 1.) * r * r
            //         + 7. / 2. * (1. - epsilon) * r
            //         + (6.125 * epsilon - 5.125);
            // };

        };
    };
    
    return wFuncValue;

};

//------------------------------------------------------------------------------
//---------------------------SETS THE WEIGHT FUNCTION---------------------------
//------------------------------------------------------------------------------
template<>
void Arlequin<2>::setWeightFunction(double val){
    
    typename Nodes::VecLocD x;
    int numIntPoints;
    double wFuncValue;

    double epsilon = 1.e-6;
    double lambda = 1.0;
 
    for (int jelCoarse = 0; jelCoarse < numElemCoarse; jelCoarse++){
        
        numIntPoints = elementsCoarse_[jelCoarse] -> 
            getNumberOfIntegrationPoints();
        
        for (int i = 0; i < numIntPoints; i++){
            
            double r = elementsCoarse_[jelCoarse]->getIntegPointDistFunction(i);

            if (r < 0){
                wFuncValue = epsilon;
            } else {
                if (r >= lambda){
                    wFuncValue = 1.;
                } else {
                    wFuncValue = (1. - epsilon) / lambda * r;
                    //wFuncValue = epsilon;
                    // wFuncValue = 1. - 3.*(1.-epsilon)/(lambda*lambda) * r * r
                    //     - 2.*(1.-epsilon)/(lambda*lambda*lambda) * r * r * r;
                    
                    if (wFuncValue < epsilon) wFuncValue = epsilon;
                };
                
            };            

            //wFuncValue = weightFunctionCoarseValue(r,epsilon);

            elementsCoarse_[jelCoarse] -> 
                setIntegPointWeightFunction(i,wFuncValue);
        };    
        
        typename Elements::Connectivity connec;

        connec = elementsCoarse_[jelCoarse] -> getConnectivity();

        for (int i = 0; i < 6; i++){        
            double r = nodesCoarse_[connec(i)] -> getDistFunction();

            if (r < 0){
                wFuncValue = epsilon;
            } else {
                if (r >= lambda){
                    wFuncValue = 1.;
                } else {
                    wFuncValue = (1. - epsilon) / lambda * r;
                    //wFuncValue = epsilon;
                    // wFuncValue = 1. - 3.*(1.-epsilon)/(lambda*lambda) * r * r
                    //     - 2.*(1.-epsilon)/(lambda*lambda*lambda) * r * r * r;

                    if (wFuncValue < epsilon) wFuncValue = epsilon;
                };
            };  
      
            //wFuncValue = weightFunctionCoarseValue(r,epsilon);
            
            elementsCoarse_[jelCoarse] -> 
                setIntegPointWeightFunction(i,wFuncValue);
                      
            //wFuncValue = weightFunctionCoarseValue(r,epsilon);

            nodesCoarse_[connec(i)] -> 
                setWeightFunction(wFuncValue);
        };         
    };


    for (int jelFine = 0; jelFine < numElemFine; jelFine++){      
        numIntPoints = elementsFine_[jelFine] -> getNumberOfIntegrationPoints();

        for (int i = 0; i < numIntPoints; i++){

            double r = elementsFine_[jelFine] -> getIntegPointDistFunction(i);

            if (r < 0){
                wFuncValue = 1.;
            }else{
                if (r >= lambda){
                    wFuncValue = epsilon;
                } else {
                    wFuncValue = 1-(1. - epsilon) / lambda * r;
                    //wFuncValue = 1. - epsilon;
                    // wFuncValue = 3.*(1.-epsilon)/(lambda*lambda) * r * r
                    //     + 2.*(1.-epsilon)/(lambda*lambda*lambda) * r * r * r;

                    if (wFuncValue > (1. - epsilon)) wFuncValue = 1. - epsilon;
                };  
            };

            elementsFine_[jelFine] -> setIntegPointWeightFunction(i,wFuncValue);
        };
    };

    for (int i=0; i<numNodesFine; i++){

        double r = nodesFine_[i] -> getDistFunction();
        
        if (r < 0){
            wFuncValue = 1.;
        }else{
            if (r >= lambda){
                wFuncValue = epsilon;
            } else {
                wFuncValue = 1- (1. - epsilon) / lambda * r;
                //wFuncValue = 1. - epsilon;
                // wFuncValue = 3.*(1.-epsilon)/(lambda*lambda) * r * r
                //     + 2.*(1.-epsilon)/(lambda*lambda*lambda) * r * r * r;
                
                if (wFuncValue > (1. - epsilon)) wFuncValue = 1. - epsilon;

            };  
        };
        
        // wFuncValue = weightFunctionFineValue(r,epsilon);

        nodesFine_[i] -> setWeightFunction(wFuncValue);
    };
    

     
    return;
};

//------------------------------------------------------------------------------
//----------------------------PRINT VELOCITY RESULTS----------------------------
//------------------------------------------------------------------------------
template<>
void Arlequin<2>::printVelocity(int step) {

    //PRINT COARSE MODEL RESULTS
    
    std::string result;
    std::ostringstream convert;

    convert << step+100000;
    result = convert.str();
    std::string s = "saidaVelCoarse"+result+".vtu";
    
    std::fstream output_v(s.c_str(), std::ios_base::out);

    output_v << "<?xml version=\"1.0\"?>" << std::endl
             << "<VTKFile type=\"UnstructuredGrid\">" << std::endl
             << "  <UnstructuredGrid>" << std::endl
             << "  <Piece NumberOfPoints=\"" << numNodesCoarse
             << "\"  NumberOfCells=\"" << numElemCoarse
             << "\">" << std::endl;

    //WRITE NODAL COORDINATES
    output_v << "    <Points>" << std::endl
             << "      <DataArray type=\"Float64\" "
             << "NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

    for (int i = 0; i < numNodesCoarse; i++){
        typename Nodes::VecLocD x;
        x = nodesCoarse_[i] -> getCoordinates();
        output_v << x(0) << " " << x(1) << " " << 0.0 << std::endl;        
    };
    output_v << "      </DataArray>" << std::endl
             << "    </Points>" << std::endl;
    
    //WRITE ELEMENT CONNECTIVITY
    output_v << "    <Cells>" << std::endl
             << "      <DataArray type=\"Int32\" "
             << "Name=\"connectivity\" format=\"ascii\">" << std::endl;
    
    for (int i = 0; i < numElemCoarse; i++){
        typename Elements::Connectivity connec;
        connec = elementsCoarse_[i] -> getConnectivity();
        output_v << connec(0) << " " << connec(1) << " " << connec(2) << " " \
                 << connec(3) << " " << connec(4) << " " << connec(5) << \
            std::endl;
    };
    output_v << "      </DataArray>" << std::endl;
  
    //WRITE OFFSETS IN DATA ARRAY
    output_v << "      <DataArray type=\"Int32\""
             << " Name=\"offsets\" format=\"ascii\">" << std::endl;
    
    int aux = 0;
    for (int i = 0; i < numElemCoarse; i++){
        output_v << aux + 6 << std::endl;
        aux += 6;
    };
    output_v << "      </DataArray>" << std::endl;
  
    //WRITE ELEMENT TYPES
    output_v << "      <DataArray type=\"UInt8\" Name=\"types\" "
             << "format=\"ascii\">" << std::endl;
    
    for (int i = 0; i < numElemCoarse; i++){
        output_v << 22 << std::endl;
    };

    output_v << "      </DataArray>" << std::endl
             << "    </Cells>" << std::endl;

    //WRITE NODAL RESULTS
    output_v << "    <PointData>" << std::endl;
    output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
             << "Name=\"Velocity\" format=\"ascii\">" << std::endl;

    for (int i=0; i<numNodesCoarse; i++){        
        output_v << nodesCoarse_[i] -> getVelocity(0) << " "              
                 << nodesCoarse_[i] -> getVelocity(1) << " " 
                 << 0. << std::endl;
    };
    output_v << "      </DataArray> " << std::endl;

    output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
             << "Name=\"Real Velocity\" format=\"ascii\">" << std::endl;

    for (int i=0; i<numNodesCoarse; i++){        
        output_v << nodesCoarse_[i] -> getVelocityArlequin(0) << " "
                 << nodesCoarse_[i] -> getVelocityArlequin(1) << " " 
                 << 0. << std::endl;
    };
    output_v << "      </DataArray> " << std::endl;

    output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
             << "Name=\"Lagrange Multipliers\" format=\"ascii\">" << std::endl;

    for (int i=0; i<numNodesCoarse; i++){        
        output_v << nodesCoarse_[i] -> getLagrangeMultiplier(0) << " "
                 << nodesCoarse_[i] -> getLagrangeMultiplier(1) << " " 
                 << 0. << std::endl;
    };
    output_v << "      </DataArray> " << std::endl;

    output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
             << "Name=\"Element\" format=\"ascii\">" << std::endl;

    for (int i=0; i<numNodesCoarse; i++){
        output_v << nodesCoarse_[i] -> getNodalElemCorrespondence() << std::endl;
    };
    output_v << "      </DataArray> " << std::endl;

    output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
             << "Name=\"Gradient\" format=\"ascii\">" << std::endl;

    for (int i=0; i<numNodesCoarse; i++){        
        output_v << nodesCoarse_[i] -> getGradientComponent(0) << " "
                 << nodesCoarse_[i] -> getGradientComponent(1) << " " 
                 << 0. << std::endl;
    };
    output_v << "      </DataArray> " << std::endl;

    output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
             << "Name=\"Dist Function\" format=\"ascii\">" << std::endl;

    for (int i=0; i<numNodesCoarse; i++){        
        output_v << nodesCoarse_[i] -> getDistFunction() << std::endl;
    };
    output_v << "      </DataArray> " << std::endl;

    output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
             << "Name=\"Weight Function\" format=\"ascii\">" << std::endl;

    for (int i=0; i<numNodesCoarse; i++){
        output_v << nodesCoarse_[i] -> getWeightFunction() << std::endl;
    };
    output_v << "      </DataArray> " << std::endl;

    output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
             << "Name=\"Pressure\" format=\"ascii\">" << std::endl;

    for (int i=0; i<numNodesCoarse; i++){
        output_v << nodesCoarse_[i] -> getPressure() << std::endl;
    };
    output_v << "      </DataArray> " << std::endl;
    output_v << "    </PointData>" << std::endl; 

    //WRITE ELEMENT RESULTS
    output_v << "    <CellData>" << std::endl;
    
    output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
             << "Name=\"Process\" format=\"ascii\">" << std::endl;
    for (int i=0; i<numElemCoarse; i++){
        output_v << domDecompCoarse.first[i] << std::endl;
    };
    output_v << "      </DataArray> " << std::endl;

    output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
             << "Name=\"Weight Function\" format=\"ascii\">" << std::endl;
    for (int i=0; i<numElemCoarse; i++){
        output_v << elementsCoarse_[i] -> getIntegPointWeightFunction(0) << std::endl;
    };
    output_v << "      </DataArray> " << std::endl;
    
    output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
             << "Name=\"Free Zone\" format=\"ascii\">" << std::endl;
     int cont=0;
    for (int i=0; i<numElemCoarse; i++){
        if (elementsFreeZone_[cont] == i){
            output_v << 1.0 << std::endl;
            cont++;
        }else{
            output_v << 0.0 << std::endl;
        };
    };
    output_v << "      </DataArray> " << std::endl;

    // output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
    //          << "Name=\"Glue Zone\" format=\"ascii\">" << std::endl;
    // cont = 0;
    // for (int i=0; i<numElemCoarse; i++){
    //     if (elementsGlueZoneCoarse_[cont] == i){
    //         output_v << 1.0 << std::endl;
    //         cont++;
    //     }else{
    //         output_v << 0.0 << std::endl;
    //     };
    // };
    // output_v << "      </DataArray> " << std::endl;

    output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
             << "Name=\"Jacobian\" format=\"ascii\">" << std::endl;

    for (int i=0; i<numElemCoarse; i++){
        output_v << elementsCoarse_[i] -> getJacobian() << std::endl;
    };
    output_v << "      </DataArray> " << std::endl;

    output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
             << "Name=\"Index\" format=\"ascii\">" << std::endl;

    for (int i=0; i<numElemCoarse; i++){
        output_v << i << std::endl;
    };
    output_v << "      </DataArray> " << std::endl;

    output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
             << "Name=\"Compressibility\" format=\"ascii\">" << std::endl;

    for (int i=0; i<numElemCoarse; i++){
        output_v << elementsCoarse_[i] -> getCompressibility() << std::endl;
    };
    output_v << "      </DataArray> " << std::endl;

    output_v << "    </CellData>" << std::endl; 

    //FINALIZE OUTPUT FILE
    output_v << "  </Piece>" << std::endl
           << "  </UnstructuredGrid>" << std::endl
           << "</VTKFile>" << std::endl;



    //PRINT FINE MODEL RESULTS
    
    std::string f = "saidaVelFine"+result+".vtu";
    
    std::fstream output_vf(f.c_str(), std::ios_base::out);

    output_vf << "<?xml version=\"1.0\"?>" << std::endl
             << "<VTKFile type=\"UnstructuredGrid\">" << std::endl
             << "  <UnstructuredGrid>" << std::endl
             << "  <Piece NumberOfPoints=\"" << numNodesFine
             << "\"  NumberOfCells=\"" << numElemFine
             << "\">" << std::endl;

    //WRITE NODAL COORDINATES
    output_vf << "    <Points>" << std::endl
             << "      <DataArray type=\"Float64\" "
             << "NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

    for (int i = 0; i < numNodesFine; i++){
        typename Nodes::VecLocD x;
        x = nodesFine_[i] -> getCoordinates();
        output_vf << x(0) << " " << x(1) << " " << 0.0 << std::endl;        
    };
    output_vf << "      </DataArray>" << std::endl
             << "    </Points>" << std::endl;
    
    //WRITE ELEMENT CONNECTIVITY
    output_vf << "    <Cells>" << std::endl
             << "      <DataArray type=\"Int32\" "
             << "Name=\"connectivity\" format=\"ascii\">" << std::endl;
    
    for (int i = 0; i < numElemFine; i++){
        typename Elements::Connectivity connec;
        connec = elementsFine_[i] -> getConnectivity();
        output_vf << connec(0) << " " << connec(1) << " " << connec(2) << " " \
                 << connec(3) << " " << connec(4) << " " << connec(5) << \
            std::endl;
    };
    output_vf << "      </DataArray>" << std::endl;
  
    //WRITE OFFSETS IN DATA ARRAY
    output_vf << "      <DataArray type=\"Int32\""
             << " Name=\"offsets\" format=\"ascii\">" << std::endl;
    
    aux = 0;
    for (int i = 0; i < numElemFine; i++){
        output_vf << aux + 6 << std::endl;
        aux += 6;
    };
    output_vf << "      </DataArray>" << std::endl;
  
    //WRITE ELEMENT TYPES
    output_vf << "      <DataArray type=\"UInt8\" Name=\"types\" "
             << "format=\"ascii\">" << std::endl;
    
    for (int i = 0; i < numElemFine; i++){
        output_vf << 22 << std::endl;
    };

    output_vf << "      </DataArray>" << std::endl
             << "    </Cells>" << std::endl;

    //WRITE NODAL RESULTS
    output_vf << "    <PointData>" << std::endl;
    output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
             << "Name=\"Velocity\" format=\"ascii\">" << std::endl;

    for (int i=0; i<numNodesFine; i++){
        output_vf << nodesFine_[i] -> getVelocity(0) << " "              
                  << nodesFine_[i] -> getVelocity(1) << " " 
                  << 0. << std::endl;
    };
    output_vf << "      </DataArray> " << std::endl;

    output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
             << "Name=\"Real Velocity\" format=\"ascii\">" << std::endl;

    for (int i=0; i<numNodesFine; i++){
        output_vf << nodesFine_[i] -> getVelocityArlequin(0) << " " 
                  << nodesFine_[i] -> getVelocityArlequin(1) << " " 
                  << 0. << std::endl;
    };
    output_vf << "      </DataArray> " << std::endl;

    output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
             << "Name=\"Lag Multipliers\" format=\"ascii\">" << std::endl;

    for (int i=0; i<numNodesFine; i++){
        output_vf << nodesFine_[i] -> getLagrangeMultiplier(0) << " " 
                  << nodesFine_[i] -> getLagrangeMultiplier(1) << " " 
                  << 0. << std::endl;
    };
    output_vf << "      </DataArray> " << std::endl;

    output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
             << "Name=\"Gradient\" format=\"ascii\">" << std::endl;

    for (int i=0; i<numNodesFine; i++){
        output_vf << nodesFine_[i] -> getGradientComponent(0) << " " 
                  << nodesFine_[i] -> getGradientComponent(1) << " " 
                  << 0. << std::endl;
    };
    output_vf << "      </DataArray> " << std::endl;

    output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
             << "Name=\"Dist Function\" format=\"ascii\">" << std::endl;

    for (int i=0; i<numNodesFine; i++){
        output_vf << nodesFine_[i] -> getDistFunction() << std::endl;
    };
    output_vf << "      </DataArray> " << std::endl;

    output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
             << "Name=\"Element\" format=\"ascii\">" << std::endl;

    for (int i=0; i<numNodesFine; i++){
        output_vf << nodesFine_[i] -> getNodalElemCorrespondence() << std::endl;
    };
    output_vf << "      </DataArray> " << std::endl;

    output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
             << "Name=\"Weight Function\" format=\"ascii\">" << std::endl;
    for (int i=0; i<numNodesFine; i++){
        output_vf << nodesFine_[i] -> getWeightFunction() << std::endl;
    };
    output_vf << "      </DataArray> " << std::endl;

    output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
             << "Name=\"Pressure\" format=\"ascii\">" << std::endl;
    for (int i=0; i<numNodesFine; i++){
        output_vf << nodesFine_[i] -> getPressure() << std::endl;
    };
    output_vf << "      </DataArray> " << std::endl;

    output_vf << "    </PointData>" << std::endl; 

    //WRITE ELEMENT RESULTS
    output_vf << "    <CellData>" << std::endl;
    
    output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
             << "Name=\"Process\" format=\"ascii\">" << std::endl;
    for (int i=0; i<numElemFine; i++){
        output_vf << domDecompFine.first[i] << std::endl;
    };
    output_vf << "      </DataArray> " << std::endl;

    output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
             << "Name=\"Glue Zone\" format=\"ascii\">" << std::endl;
    cont=0;
    for (int i=0; i<numElemFine; i++){
        if (elementsGlueZoneFine_[cont] == i){
            output_vf << 1.0 << std::endl;
            cont += 1; 
        }else{
            output_vf << 0.0 << std::endl;
        };
    };
    output_vf << "      </DataArray> " << std::endl;
    
    output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
             << "Name=\"Weight Function\" format=\"ascii\">" << std::endl;
    for (int i=0; i<numElemFine; i++){
        output_vf << elementsFine_[i] -> getIntegPointWeightFunction(0) << std::endl;
    };
    output_vf << "      </DataArray> " << std::endl;

    output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
             << "Name=\"Jacobian\" format=\"ascii\">" << std::endl;

    for (int i=0; i<numElemFine; i++){
        output_vf << elementsFine_[i] -> getJacobian() << std::endl;
    };
    output_vf << "      </DataArray> " << std::endl;

    output_vf << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
             << "Name=\"Compressibility\" format=\"ascii\">" << std::endl;

    for (int i=0; i<numElemFine; i++){
        output_vf << elementsFine_[i] -> getCompressibility() << std::endl;
    };
    output_vf << "      </DataArray> " << std::endl;

    output_vf << "    </CellData>" << std::endl; 

    //FINALIZE OUTPUT FILE
    output_vf << "  </Piece>" << std::endl
           << "  </UnstructuredGrid>" << std::endl
           << "</VTKFile>" << std::endl;



    std::string c = "saidaCoupling"+result+".vtu";
    
    std::fstream output_c(c.c_str(), std::ios_base::out);

    output_c << "<?xml version=\"1.0\"?>" << std::endl
             << "<VTKFile type=\"UnstructuredGrid\">" << std::endl
             << "  <UnstructuredGrid>" << std::endl
             << "  <Piece NumberOfPoints=\"" << numNodesGlueZoneFine
             << "\"  NumberOfCells=\"" << numElemGlueZoneFine
             << "\">" << std::endl;

    //WRITE NODAL COORDINATES
    output_c << "    <Points>" << std::endl
             << "      <DataArray type=\"Float64\" "
             << "NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

    for (int i = 0; i < numNodesGlueZoneFine; i++){
        typename Nodes::VecLocD x;
        x = nodesLagrangeFine_[i] -> getCoordinates();
        output_c << x(0) << " " << x(1) << " " << 0.0 << std::endl;        
    };
    output_c << "      </DataArray>" << std::endl
             << "    </Points>" << std::endl;
    
    //WRITE ELEMENT CONNECTIVITY
    output_c << "    <Cells>" << std::endl
             << "      <DataArray type=\"Int32\" "
             << "Name=\"connectivity\" format=\"ascii\">" << std::endl;
    
    for (int i = 0; i < numElemGlueZoneFine; i++){
        typename Elements::Connectivity connec;
        connec = glueZoneFine_[i] -> getConnectivity();
        output_c << connec(0) << " " << connec(1) << " " << connec(2) << " " \
                 << connec(3) << " " << connec(4) << " " << connec(5) << \
            std::endl;
    };
    output_c << "      </DataArray>" << std::endl;
  
    //WRITE OFFSETS IN DATA ARRAY
    output_c << "      <DataArray type=\"Int32\""
             << " Name=\"offsets\" format=\"ascii\">" << std::endl;
    
    aux = 0;
    for (int i = 0; i < numElemGlueZoneFine; i++){
        output_c << aux + 6 << std::endl;
        aux += 6;
    };
    output_c << "      </DataArray>" << std::endl;
  
    //WRITE ELEMENT TYPES
    output_c << "      <DataArray type=\"UInt8\" Name=\"types\" "
             << "format=\"ascii\">" << std::endl;
    
    for (int i = 0; i < numElemGlueZoneFine; i++){
        output_c << 22 << std::endl;
    };

    output_c << "      </DataArray>" << std::endl
             << "    </Cells>" << std::endl;

    //WRITE NODAL RESULTS
    output_c << "    <PointData>" << std::endl;
    output_c << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
             << "Name=\"Velocity\" format=\"ascii\">" << std::endl;

    for (int i=0; i<numNodesGlueZoneFine; i++){
        output_c << nodesFine_[i] -> getVelocity(0) << " "              
                  << nodesFine_[i] -> getVelocity(1) << " " 
                  << 0. << std::endl;
    };
    output_c << "      </DataArray> " << std::endl;


    output_c << "    </PointData>" << std::endl; 

    //WRITE ELEMENT RESULTS
    output_c << "    <CellData>" << std::endl;
    

    output_c << "    </CellData>" << std::endl; 

    //FINALIZE OUTPUT FILE
    output_c << "  </Piece>" << std::endl
           << "  </UnstructuredGrid>" << std::endl
           << "</VTKFile>" << std::endl;
};


//------------------------------------------------------------------------------
//------------SETS COARSE/FINE MESHES AND GETS ITS BASIC INFORMATIONS-----------
//------------------------------------------------------------------------------
template<>
void Arlequin<2>::setFluidModels(FluidMesh coarse, FluidMesh fine){

    coarseModel = coarse;
    fineModel = fine;

    //Gets Fine and Coarse models basic information
    numElemCoarse = coarseModel.elements_.size();
    numElemFine   = fineModel.elements_.size();
    numNodesCoarse = coarseModel.nodes_.size();
    numNodesFine   = fineModel.nodes_.size();

    nodesCoarse_  = coarseModel.nodes_;
    nodesFine_    = fineModel.nodes_;

    elementsCoarse_ = coarseModel.elements_;
    elementsFine_   = fineModel.elements_;

    boundaryCoarse_ = coarseModel.boundary_;
    boundaryFine_   = fineModel.boundary_;

    domDecompCoarse = coarseModel.getDomainDecomposition();
    domDecompFine = fineModel.getDomainDecomposition();

    numTimeSteps = fineModel.getNumberOfTimeSteps();
    dTime = fineModel.getTimeStep();

    for (int i=0; i < numElemFine; i++){
        elementsFine_[i] -> setModel(true);
    };
    for (int i=0; i < numElemCoarse; i++){
        elementsCoarse_[i] -> setModel(false);
    };

    //Sets the element boxes for all elements in both coarse and fine models
    setElementBoxes();
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);      

};

//------------------------------------------------------------------------------
//----------------COMPUTE ARLEQUIN COUPLED MOVING LAPLACE PROBLEM---------------
//------------------------------------------------------------------------------
template<>
int Arlequin<2>::solveSteadyArlequinMovingLaplaceProblem(int iterNumber, 
                                                         double tolerance,
                                                         int steps){

    Mat               A;
    Vec               b, u, All;
    PetscErrorCode    ierr;
    PetscInt          Istart, Iend, Ii, Ione, iterations;
    KSP               ksp;
    PC                pc;
    VecScatter        ctx;
    PetscScalar       val;
    
    //Construct the glue zone based on some defined criterion
    setCouplingZone();

    //Computes the Weight function for all the finite elements
    setWeightFunction(4.);

    //Computes the Nodal correspondence between fine nodes and coarse elements
    setNodalCorrespondenceFine();

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);      

    // Computes the system size
    int sysSize = numNodesCoarse + numNodesFine + numNodesFine;

    for (int iSteps = 0; iSteps < steps; iSteps++){


        for (int i=0; i<numNodesFine; i++){
            
            typename Nodes::VecLocD x;

            x = nodesFine_[i] -> getCoordinates();

            x(0) += 0.1;
            //            std::cout << "x0 " << x(0) << std::endl;
            nodesFine_[i] -> setCoordinates(x);
        };


        
        for (int i=0; i<numNodesFine; i++){
            nodesFine_[i] -> clearVariables();
        };
        for (int i=0; i<numNodesCoarse; i++){
            nodesCoarse_[i] -> clearVariables();
        };
        for (int i=0; i<numElemFine; i++){
            elementsFine_[i] -> clearVariables();
        };
        for (int i=0; i<numElemCoarse; i++){
            elementsCoarse_[i] -> clearVariables();
        }; 
                
        //Sets the element boxes for all elements in both coarse and fine models
        setElementBoxes();
        
        //Computes the Nodal correspondence between fine nodes and coarse elements
        setNodalCorrespondenceFine();

        // //Construct the glue zone based on some defined criterion
        // setCouplingZoneFine();

        // //Computes the Weight function for all the finite elements
        setWeightFunction(4. + 0.1 * (iSteps+1));


    for (int inewton = 0; inewton < iterNumber; inewton++){
        
        // Preallocates the matrix
        ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                            sysSize, sysSize, 500, NULL, 500, NULL, &A); 
        CHKERRQ(ierr);
        
        // Divides the matrix between the processes
        ierr = MatGetOwnershipRange(A, &Istart, &Iend);CHKERRQ(ierr);
        
        //Create PETSc vectors
        ierr = VecCreate(PETSC_COMM_WORLD, &b); CHKERRQ(ierr);
        ierr = VecSetSizes(b, PETSC_DECIDE, sysSize); CHKERRQ(ierr);
        ierr = VecSetFromOptions(b); CHKERRQ(ierr); 
        ierr = VecDuplicate(b, &u); CHKERRQ(ierr);
        ierr = VecDuplicate(b, &All); CHKERRQ(ierr);
        
        //std::cout << "Istart = " << Istart << " Iend = " << Iend << std::endl;
        
        //----------------------------------------------------------------------
        //--------------------BEGIN OF LINEAR SYSTEM ASSEMBLY-------------------
        //----------------------------------------------------------------------

        //Coarse mesh
        for (int jel = 0; jel < numElemCoarse; jel++){   
            
            if (domDecompCoarse.first[jel] == rank) {
                
                //Compute Element matrix
                elementsCoarse_[jel] -> getSteadyLaplace();
                                
                typename Elements::LocalMatrix Ajac;
                typename Elements::LocalVector Rhs;
                typename Elements::Connectivity connec;
                
                //Gets element connectivity, jacobian and rhs 
                connec = elementsCoarse_[jel] -> getConnectivity();
                Ajac = elementsCoarse_[jel] -> getJacNRMatrix();
                Rhs = elementsCoarse_[jel] -> getRhsVector();
                
                //Disperse local contributions into the global matrix
                for (int i=0; i<6; i++){
                    for (int j=0; j<6; j++){
                        if (fabs(Ajac(2*i  ,2*j  )) >= 1.e-8){
                            int dof_i = connec(i);
                            int dof_j = connec(j);
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    
                                                &Ajac(2*i  ,2*j  ),ADD_VALUES);
                        };
                    };                                       
                    //Rhs vector
                    if (fabs(Rhs(2*i  )) >= 1.e-8){
                        int dof_i = connec(i);
                        ierr = VecSetValues(b,1,&dof_i,&Rhs(2*i  ),ADD_VALUES);
                    };
                };
            };
        };
        
        //Fine mesh
        for (int jel = 0; jel < numElemFine; jel++){   
            
            if (domDecompFine.first[jel] == rank) {
                
                //Compute Element matrix
                elementsFine_[jel] -> getSteadyLaplace();
                
                typename Elements::LocalMatrix Ajac;
                typename Elements::LocalVector Rhs;
                typename Elements::Connectivity connec;
                
                //Gets element connectivity, jacobian and rhs 
                connec = elementsFine_[jel] -> getConnectivity();
                Ajac = elementsFine_[jel] -> getJacNRMatrix();
                Rhs = elementsFine_[jel] -> getRhsVector();
                
                //Disperse local contributions into the global matrix
                for (int i=0; i<6; i++){
                    for (int j=0; j<6; j++){
                        if (fabs(Ajac(2*i  ,2*j  )) >= 1.e-8){
                            int dof_i = numNodesCoarse + connec(i);
                            int dof_j = numNodesCoarse + connec(j);
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    \
                                                &Ajac(2*i  ,2*j  ),ADD_VALUES);
                        };
                    };
                                        
                    //Rhs vector
                    if (fabs(Rhs(2*i  )) >= 1.e-8){
                        int dof_i = numNodesCoarse + connec(i);
                        ierr = VecSetValues(b,1,&dof_i,&Rhs(2*i  ),ADD_VALUES);
                    };
                };
            };
        };
                
        //Lagrange Multipliers
        for (int l=0; l< numElemGlueZoneFine; l++){
          
            int jel = elementsGlueZoneFine_[l];
            
            if (domDecompFine.first[jel] == rank) {
                
                typename Elements::LocalMatrix Ajac;
                typename Elements::LocalVector Rhs;
                typename Elements::Connectivity connec, connecC;
                
                connec = elementsFine_[jel] -> getConnectivity();
                
                // FINE MESH

                //Computes element matrix
                elementsFine_[jel] -> getLagrangeMultipliersSameMesh();
                
                //Gets element matrice and rhs vectors
                Ajac = - elementsFine_[jel] -> getJacNRMatrix();
                Rhs = - elementsFine_[jel] -> getRhsVector();
                
                typename Elements::LocalVector rhsLagMult,U_;
                
                U_.clear();
                for (int i = 0; i < 6; i++){
                    U_(2*i) = nodesFine_[connec(i)] -> getLagrangeMultiplier(0);
                };
                
                noalias(rhsLagMult) = - prod(trans(Ajac),U_);
                
                //Disperse local contributions into the global matrix
                for (int i=0; i<6; i++){
                    for (int j=0; j<6; j++){
                        if (fabs(Ajac(2*i  ,2*j  )) >= 1.e-8){
                            int d_i = numNodesCoarse + numNodesFine + connec(i);
                            int d_j = numNodesCoarse + connec(j);
                                    
                            ierr = MatSetValues(A,1,&d_i,1,&d_j,
                                                &Ajac(2*i  ,2*j  ),ADD_VALUES);
                            ierr = MatSetValues(A,1,&d_j,1,&d_i,
                                                &Ajac(2*j  ,2*i  ),ADD_VALUES);
                        };
                    };
                    //Rhs vector
                    if (fabs(Rhs(2*i  )) >= 1.e-8){
                        int dof_i = numNodesCoarse + numNodesFine + connec(i);
                        ierr = VecSetValues(b,1,&dof_i,&Rhs(2*i  ), ADD_VALUES);
                    };
                    if (fabs(rhsLagMult(2*i  )) >= 1.e-8){
                        int dof_i = numNodesCoarse + connec(i);
                        ierr = VecSetValues(b,1,&dof_i,&rhsLagMult(2*i  )
                                            ,ADD_VALUES);
                    };
                };      
                
                //COAESE MESH
                
                //Counts number of coarse mesh intersecting the 
                //fine element
                int numberIntPoints = elementsFine_[jel] -> 
                    getNumberOfIntegrationPoints();
                int aux;
                
                std::vector<int> ele, diffElem;
                ele.clear();
                diffElem.clear();
                for (int i=0; i<numberIntPoints; i++){
                    aux = elementsFine_[jel] -> 
                        getIntegPointCorrespondenceElement(i);
                    ele.push_back(aux);
                };
                
                int numElemIntersect = 1;
                int flag = 0;
                diffElem.push_back(ele[0]);
                
                for (int i = 1; i<numberIntPoints; i++){
                    flag = 0;
                    for (int j = 0; j<numElemIntersect; j++){
                        if (ele[i] == diffElem[j]) {
                            break;
                        }else{
                            flag++;
                        };
                        if(flag == numElemIntersect){
                            numElemIntersect++;
                            diffElem.push_back(ele[i]);
                        };
                    };
                };
                
                //Compute the Lagrange Multiplier element matrix
                for (int ielem = 0; ielem < numElemIntersect; ielem++){
                    
                    int iElemCoarse = diffElem[ielem];
                    
                    elementsFine_[jel] -> 
                        getLagrangeMultipliersDifferentMesh(iElemCoarse);
                    
                    //Computes element matrix
                    Ajac = elementsFine_[jel] -> getJacNRMatrix();
                    
                    typename Elements::LocalVector rhsLagMult,U_;
                    
                    U_.clear();
                    for (int i = 0; i < 6; i++){
                        U_(2*i  ) = 
                            nodesFine_[connec(i)] -> getLagrangeMultiplier(0);
                    };
                    
                    noalias(rhsLagMult) = -prod(trans(Ajac),U_);
                    
                    std::pair<Elements::LocalVector, 
                              Elements::LocalMatrix>  lagMult;
                    
                    connecC = elementsCoarse_[iElemCoarse] -> getConnectivity();
                            
                    lagMult = elementsCoarse_[iElemCoarse] -> 
                        getRhsVectorAndBoundaryConditions(Ajac);
                            
                    //Disperse local contribution into the global matrix
                    for (int i=0; i<6; i++){
                        for (int j=0; j<6; j++){
                            if (fabs(lagMult.second(2*i,2*j)) >= 1.e-8){
                                int dof_i = numNodesCoarse + numNodesFine 
                                    + connec(i);
                                int dof_j = connecC (j);
                                
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,
                                                    &lagMult.second(2*i  ,2*j  )
                                                    ,ADD_VALUES);
                                
                                ierr = MatSetValues(A,1,&dof_j,1,&dof_i,
                                                    &lagMult.second(2*i  ,2*j  )
                                                    ,ADD_VALUES);
                            };
                        };
                        //Rhs vector
                        if (fabs(lagMult.first(2*i  )) >= 1.e-8){
                            int d_i = numNodesCoarse +numNodesFine + connec(i);
                            ierr = VecSetValues(b,1,&d_i,&lagMult.first(2*i  ),
                                                ADD_VALUES);
                        };
                        
                        if (fabs(rhsLagMult(2*i  )) >= 1.e-8){
                            int dof_i = connecC(i);
                            ierr = VecSetValues(b,1,&dof_i,&rhsLagMult(2*i  ),
                                                ADD_VALUES);
                        };
                    };                                 
                }; //Number of intersections
            }; // if element belongs to the glue zone
        }; // Glue zone
        //----------------------------------------------------------------------
        //---------------------END OF LINEAR SYSTEM ASSEMBLY--------------------
        //----------------------------------------------------------------------
        
        //Assemble matrices and vectors
        ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        
        ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
        ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
        
        //ierr = MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
        //ierr = VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
        
        //Create KSP context to solve the linear system
        ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
        
        ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
        
        ierr = KSPSetTolerances(ksp,1.e-7,1.e-10,PETSC_DEFAULT,
                                10000);CHKERRQ(ierr);
        
        ierr = KSPGMRESSetRestart(ksp, 10); CHKERRQ(ierr);
      
        ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
        
        ierr = PCSetType(pc, PCJACOBI);CHKERRQ(ierr);
        
        //ierr = KSPSetPCSide(ksp, PC_RIGHT);
        //ierr = KSPSetType(ksp,KSPTFQMR); CHKERRQ(ierr);

        ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
        //ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
        
        ierr = KSPSolve(ksp,b,u);CHKERRQ(ierr);
        
        ierr = KSPGetTotalIterations(ksp, &iterations);

        std::cout << "GMRES Iterations = " << iterations << std::endl;
        
        //ierr = VecView(u,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
        
        //Gathers the solution vector to the master process
        ierr = VecScatterCreateToAll(u, &ctx, &All);CHKERRQ(ierr);
        
        ierr = VecScatterBegin(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);
        CHKERRQ(ierr);
        
        ierr = VecScatterEnd(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);
        CHKERRQ(ierr);
        
        ierr = VecScatterDestroy(&ctx);CHKERRQ(ierr);
                
        //Updates nodal values
        double u_[2];
        double norm = 0.;
        Ione = 1;

        for (int i = 0; i < numNodesCoarse; ++i){
            Ii = i;
            ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
            u_[0] = val;
            norm += val*val;
            nodesCoarse_[i] -> incrementVelocity(0,u_[0]);
        };

        for (int i = 0; i < numNodesFine; ++i){
            Ii = numNodesCoarse + i;
            ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
            u_[0] = val;
            norm += val*val;
            nodesFine_[i] -> incrementVelocity(0,u_[0]);
        };

        for (int i = 0; i < numNodesFine; ++i){
            Ii = numNodesCoarse + numNodesFine + i;
            ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
            u_[0] = val;
            nodesFine_[i] -> incrementLagrangeMultiplier(0,u_[0]);
        };
        
        //Computes the solution vector norm
        ierr = VecNorm(u,NORM_2,&val);CHKERRQ(ierr);



        
        QuadShapeFunction<2>                       shapeQuad;
        typename QuadShapeFunction<2>::Values      phi_;

        for (int i = 0; i<numNodesFine; i++){
            
            typename Nodes::VecLocD xsi;
            typename Quadrature::NodalValuesQuad p_coarse;
            typename Elements::Connectivity connecCoarse;

            double val = 0.;
            
            int elCoarse = nodesFine_[i] -> getNodalElemCorrespondence();
            xsi = nodesFine_[i] -> getNodalXsiCorrespondence();
            
            connecCoarse = elementsCoarse_[elCoarse] -> getConnectivity();
            
            for (int j=0; j<6; j++){
                p_coarse(j) = nodesCoarse_[connecCoarse(j)] -> getVelocity(0);
            };

            shapeQuad.evaluate(xsi,phi_);

            for (int j=0; j<6; j++){
                val += p_coarse(j) * phi_(j);
            };
            
            double wFunc = nodesFine_[i] -> getWeightFunction();

            double pot = nodesFine_[i] -> getVelocity(0) * wFunc 
                + val * (1 - wFunc);

            nodesFine_[i] -> setPotential(pot);
            
        };

        for (int i = 0; i<numNodesCoarse; i++){
            nodesCoarse_[i] -> setPotential(nodesCoarse_[i] -> getVelocity(0));
        };
        
        if(rank == 0){
            std::cout << "STEP  = " << iSteps << std::endl;
            
            std::cout << "Iteration = " << inewton 
                      << "      Du Norm = " << val 
                      << " " << sqrt(norm) << std::endl;
        };

        ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
        ierr = VecDestroy(&b); CHKERRQ(ierr);
        ierr = VecDestroy(&u); CHKERRQ(ierr);
        ierr = VecDestroy(&All); CHKERRQ(ierr);
        ierr = MatDestroy(&A); CHKERRQ(ierr);


        if(val <= tolerance){
            break;            
        };             
    };




    //Postprocessing results

    for (int i=0; i<numElemFine; i++){
        elementsFine_[i] -> computeNodalGradient();
    };
    
    for (int i=0; i<numElemCoarse; i++){
        elementsCoarse_[i] -> computeNodalGradient();
    };        
    


    // if (rank == 0){
    //     int numPoints = 100;
    //     resultsX_.reserve(numPoints); resultsY_.reserve(numPoints);
        
    //     typename Nodes::VecLocD x;
        
    //     for (int i=0; i<numPoints; i++){
    //         x(0) = 4.0;
    //         x(1) = i / 33. + 4.5;
    //         Nodes *node = new Nodes(x, i);
    //         resultsY_.push_back(node);
    //         x(0) = i / 33. + 0.5;
    //         x(1) = 4.;
    //         Nodes *node2 = new Nodes(x, i);
    //         resultsX_.push_back(node2);
    //     };

    //     std::pair<int, ublas::bounded_vector<double,2> > corresp;

    //     for (int i=0; i<numPoints; i++){
    //         x = resultsX_[i] -> getCoordinates();
            
    //         corresp = searchNodeCorrespondence(x, nodesFine_, elementsFine_, 
    //                                            elementsFine_.size());
    //         resultsX_[i] -> setNodalCorrespondence(corresp.first, 
    //                                                corresp.second);

    //         x = resultsY_[i] -> getCoordinates();
            
    //         corresp = searchNodeCorrespondence(x, nodesFine_, elementsFine_, 
    //                                            elementsFine_.size());
    //         resultsY_[i] -> setNodalCorrespondence(corresp.first, 
    //                                                corresp.second);
    //     };
        
    //     std::pair<double, ublas::bounded_vector<double,2> > resul;
    //     int el;
    //     double u[2];

    //     for (int i=0; i<numPoints; i++){

    //         x = resultsX_[i] -> getNodalXsiCorrespondence();
    //         el = resultsX_[i] -> getNodalElemCorrespondence();
            
    //         resul = elementsFine_[el] -> getPotentialResultsInPoint(x);
            
    //         resultsX_[i] -> setPotential(resul.first);
    //         u[0] = resul.second(0);
    //         u[1] = resul.second(1);
    //         resultsX_[i] -> setVelocity(u);

    //         x = resultsY_[i] -> getNodalXsiCorrespondence();
    //         el = resultsY_[i] -> getNodalElemCorrespondence();
            
    //         resul = elementsFine_[el] -> getPotentialResultsInPoint(x);
            
    //         resultsY_[i] -> setPotential(resul.first);
    //         u[0] = resul.second(0);
    //         u[1] = resul.second(1);
    //         resultsY_[i] -> setVelocity(u);
    //     };

    //     std::string ss = "resultsCubic0.5.dat";
        
    //     std::fstream out(ss.c_str(), std::ios_base::out);

    //     for (int i=0; i<numPoints; i++){
    //         out << std::scientific << resultsX_[i] -> getPotential() << " "  
    //             << resultsX_[i] -> getVelocity(0) << " "  
    //             << resultsX_[i] -> getVelocity(1) << " "  
    //             << resultsY_[i] -> getPotential() << " "  
    //             << resultsY_[i] -> getVelocity(0) << " "  
    //             << resultsY_[i] -> getVelocity(1) << std::endl;        
            
    //     };
    // };

    if (rank == 0) {
        //Computing velocity divergent
        printVelocity(iSteps);
    };

    };
    
    std::cout << "Degrees of freedom = "
              << numNodesCoarse + numNodesFine + numNodesGlueZoneFine 
              << std::endl;

    return 0;

};


//------------------------------------------------------------------------------
//----------------COMPUTE ARLEQUIN COUPLED NAVIER-STOKES PROBLEM----------------
//------------------------------------------------------------------------------
template<>
int Arlequin<2>::solveArlequinProblem(int iterNumber, double tolerance,
                                      int problem_type, int time_dependency){

    Mat               A;
    Vec               b, u, All;
    PetscErrorCode    ierr;
    PetscInt          Istart, Iend, Ii, Ione, iterations;
    KSP               ksp; 
    PC                pc;
    VecScatter        ctx;
    PetscScalar       val;
    //MatNullSpace      nullsp;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);      

    if ((problem_type < 1) || (problem_type > 2)){
        std::cout << "WRONG PROBLEM TYPE." << std::endl;
        return 0;
    };

    if (time_dependency == 0) numTimeSteps = 1;
    
    iTimeStep = 0.;

    //Construct the glue zone based on some defined criterion
    setCouplingZone();
    
    std::cout << "AQUI saiu " << rank << std::endl;     
    MPI_Barrier(PETSC_COMM_WORLD);      

    //Computes the Weight function for all the finite elements
    setWeightFunction(4.);

    //Computes the Nodal correspondence between fine nodes and coarse elements
    setNodalCorrespondenceFine();

    //Construct the free zone in coarse model
    setFreeZone();

    if(rank == 0) printVelocity(0);
    MPI_Barrier(PETSC_COMM_WORLD);      


    for (int jel = 0; jel < numElemGlueZoneFine; jel++){
        int iel = elementsGlueZoneFine_[jel];
        //elementsFine_[iel] -> setCompressibility();
    };

    // Computes the system size
    int sysSize = 3 * numNodesCoarse +  
        + 3 * numNodesFine + 2 * numNodesGlueZoneFine;
    
    for (iTimeStep = 0; iTimeStep < numTimeSteps; iTimeStep++){
        
        if (rank == 0) {std::cout << "------------------------- TIME STEP = "
                                  << iTimeStep << " -------------------------"
                                  << std::endl;}
        
        //Updates velocity and acceleration
        for (int i = 0; i < numNodesCoarse; i++){
            double accel[2], u[2], uprev[2];
            
            //Compute acceleration
            u[0] = nodesCoarse_[i] -> getVelocity(0);
            u[1] = nodesCoarse_[i] -> getVelocity(1);
            
            uprev[0] = nodesCoarse_[i] -> getPreviousVelocity(0);
            uprev[1] = nodesCoarse_[i] -> getPreviousVelocity(1);
            
            accel[0] = (u[0] - uprev[0]) / dTime;
            accel[1] = (u[1] - uprev[1]) / dTime;
            
            nodesCoarse_[i] -> setAcceleration(accel);
            
            //Updates velocity
            nodesCoarse_[i] -> setPreviousVelocity(u);
            
            //Computes predictor
            // u[0] += dTime * accel[0];
            // u[1] += dTime * accel[1];
            
            // nodesCoarse_[i] -> setVelocity(u);
            
            // //Updates acceleration
            // nodesCoarse_[i] -> setPreviousAcceleration(accel);
            
        };
        for (int i = 0; i < numNodesFine; i++){
            double accel[2], u[2], uprev[2];
            
            //Compute acceleration
            u[0] = nodesFine_[i] -> getVelocity(0);
            u[1] = nodesFine_[i] -> getVelocity(1);
            
            uprev[0] = nodesFine_[i] -> getPreviousVelocity(0);
            uprev[1] = nodesFine_[i] -> getPreviousVelocity(1);
            
            accel[0] = (u[0] - uprev[0]) / dTime;
            accel[1] = (u[1] - uprev[1]) / dTime;
            
            nodesFine_[i] -> setAcceleration(accel);
            
            //Updates velocity
            nodesFine_[i] -> setPreviousVelocity(u);
            
            //Computes predictor
            // u[0] += dTime * accel[0];
            // u[1] += dTime * accel[1];
            
            // nodesFine_[i] -> setVelocity(u);
            
            // //Updates acceleration
            // nodesFine_[i] -> setPreviousAcceleration(accel);
            
        };

        //STARTS NEWTON-RAPHSON
        for (int inewton = 0; inewton < iterNumber; inewton++){
            
            // if (iTimeStep < 100){
            //     for (int i = 0; i < numElemFine; i++){
            //         double visc = 0.0002;//*iTimeStep/100.;
            //         elementsFine_[i] -> setViscosity(visc);
            //     };
            //     for (int i = 0; i < numElemCoarse; i++){
            //         double visc = 0.0002;//*iTimeStep/100.;
            //         elementsCoarse_[i] -> setViscosity(visc);
            //     };
            // }else{
            //     for (int i = 0; i < numElemFine; i++){
            //         double visc = 0.0002;
            //         elementsFine_[i] -> setViscosity(visc);
            //     };
            //     for (int i = 0; i < numElemCoarse; i++){
            //         double visc = 0.0002;
            //         elementsCoarse_[i] -> setViscosity(visc);
            //     };
            // };


            boost::posix_time::ptime t1 =
                boost::posix_time::microsec_clock::local_time();
            
            // Preallocates the matrix
            ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                                sysSize, sysSize, 120, NULL, 600, NULL, &A); 
            CHKERRQ(ierr);
            
            // Divides the matrix between the processes
            ierr = MatGetOwnershipRange(A, &Istart, &Iend);CHKERRQ(ierr);
            
            //Create PETSc vectors
            ierr = VecCreate(PETSC_COMM_WORLD, &b); CHKERRQ(ierr);
            ierr = VecSetSizes(b, PETSC_DECIDE, sysSize); CHKERRQ(ierr);
            ierr = VecSetFromOptions(b); CHKERRQ(ierr); 
            ierr = VecDuplicate(b, &u); CHKERRQ(ierr);
            ierr = VecDuplicate(b, &All); CHKERRQ(ierr);
                        
            //------------------------------------------------------------------
            //------------------BEGIN OF LINEAR SYSTEM ASSEMBLY-----------------
            //------------------------------------------------------------------

            //Coarse mesh
            
            for (int jel = 0; jel < numElemCoarse; jel++){   
                
                if (domDecompCoarse.first[jel] == rank) {
                    
                    //Compute Element matrix
                    if (time_dependency == 0){ 
                        if (problem_type == 1){
                            elementsCoarse_[jel] -> getSteadyStokes();
                        }else{
                            if (inewton == 0){
                                elementsCoarse_[jel] -> getSteadyStokes();
                            }else{
                                elementsCoarse_[jel] -> getSteadyNavierStokes();
                            };
                        };
                    }else{
                        if (problem_type == 1){
                            elementsCoarse_[jel] -> getTransientNavierStokes();
                        }else{
                            elementsCoarse_[jel] -> getTransientNavierStokes();
                        };
                    };

                    typename Elements::LocalMatrix Ajac;
                    typename Elements::LocalVector Rhs;
                    typename Elements::Connectivity connec;
                    
                    //Gets element connectivity, jacobian and rhs 
                    connec = elementsCoarse_[jel] -> getConnectivity();
                    Ajac = elementsCoarse_[jel] -> getJacNRMatrix();
                    Rhs = elementsCoarse_[jel] -> getRhsVector();
                    
                    //Disperse local contributions into the global matrix
                    for (int i=0; i<6; i++){
                        for (int j=0; j<6; j++){
                            
                            if (fabs(Ajac(2*i  ,2*j  )) >= 1.e-15){
                                int dof_i = 2 * connec(i);
                                int dof_j = 2 * connec(j);
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    
                                                    &Ajac(2*i  ,2*j  ),
                                                    ADD_VALUES);
                            };
                            if (fabs(Ajac(2*i+1,2*j  )) >= 1.e-15){
                                int dof_i = 2 * connec(i) + 1;
                                int dof_j = 2 * connec(j);
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    
                                                    &Ajac(2*i+1,2*j  ),
                                                    ADD_VALUES);
                            };

                            if (fabs(Ajac(2*i  ,2*j+1)) >= 1.e-15){
                                int dof_i = 2 * connec(i);
                                int dof_j = 2 * connec(j) + 1;
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    
                                                    &Ajac(2*i  ,2*j+1),
                                                    ADD_VALUES);
                            };

                            if (fabs(Ajac(2*i+1,2*j+1)) >= 1.e-15){
                                int dof_i = 2 * connec(i) + 1;
                                int dof_j = 2 * connec(j) + 1;
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    
                                                    &Ajac(2*i+1,2*j+1),
                                                    ADD_VALUES);
                            };

                            //Matrix Q and Qt
                            if (fabs(Ajac(2*i  ,12+j)) >= 1.e-15){
                                int dof_i = 2 * connec(i);
                                int dof_j = 2 * numNodesCoarse + connec(j);
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,
                                                    &Ajac(2*i  ,12+j),
                                                    ADD_VALUES);
                            };
                            
                            if (fabs(Ajac(12+j,2*i  )) >= 1.e-15){
                                int dof_i = 2 * connec(i);
                                int dof_j = 2 * numNodesCoarse + connec(j);
                                ierr = MatSetValues(A,1,&dof_j,1,&dof_i,    
                                                    &Ajac(12+j,2*i  ),
                                                    ADD_VALUES);
                            };

                            if (fabs(Ajac(2*i+1,12+j)) >= 1.e-15){
                                int dof_i = 2 * connec(i) + 1;
                                int dof_j = 2 * numNodesCoarse + connec(j);
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,
                                                    &Ajac(2*i+1,12+j),
                                                    ADD_VALUES);
                            };

                            if (fabs(Ajac(12+j,2*i+1)) >= 1.e-15){
                                int dof_i = 2 * connec(i) + 1;
                                int dof_j = 2 * numNodesCoarse + connec(j);
                                ierr = MatSetValues(A,1,&dof_j,1,&dof_i,    
                                                    &Ajac(12+j,2*i+1),
                                                    ADD_VALUES);
                            };

                            if (fabs(Ajac(12+i,12+j)) >= 1.e-15){
                                int dof_i = 2 * numNodesCoarse + connec(i);
                                int dof_j = 2 * numNodesCoarse + connec(j);
                                ierr = MatSetValues(A,1,&dof_j,1,&dof_i,    
                                                    &Ajac(12+i,12+j),
                                                    ADD_VALUES);
                            };
                        };                  
                        //Rhs vector
                        if (fabs(Rhs(2*i  )) >= 1.e-15){
                            int dof_i = 2 * connec(i);
                            ierr = VecSetValues(b,1,&dof_i,&Rhs(2*i  ),
                                                ADD_VALUES);
                        };

                        if (fabs(Rhs(2*i+1)) >= 1.e-15){
                            int dof_i = 2 * connec(i) + 1;
                            ierr = VecSetValues(b,1,&dof_i,&Rhs(2*i+1),
                                                ADD_VALUES);
                        };

                        if (fabs(Rhs(12+i)) >= 1.e-15){
                            int dof_i = 2 * numNodesCoarse + connec(i);
                            ierr = VecSetValues(b,1,&dof_i,&Rhs(12+i),
                                                ADD_VALUES);
                        };
                    };
                };
            };

            //Fine mesh
            for (int jel = 0; jel < numElemFine; jel++){   
                
                if (domDecompFine.first[jel] == rank) {
                    
                    //Compute Element matrix
                    if (time_dependency == 0) {
                        if (problem_type == 1){
                            elementsFine_[jel] -> getSteadyStokes();
                        }else{
                            if (inewton == 0){
                                elementsFine_[jel] -> getSteadyStokes();
                            }else{
                                elementsFine_[jel] -> getSteadyNavierStokes();
                            };
                        };
                    }else{
                        if (problem_type == 1){
                            elementsFine_[jel] -> getTransientNavierStokes();
                        }else{
                            elementsFine_[jel] -> getTransientNavierStokes();
                        };
                    };

                    typename Elements::LocalMatrix Ajac;
                    typename Elements::LocalVector Rhs;
                    typename Elements::Connectivity connec;
                    
                    //Gets element connectivity, jacobian and rhs 
                    connec = elementsFine_[jel] -> getConnectivity();
                    Ajac = elementsFine_[jel] -> getJacNRMatrix();
                    Rhs = elementsFine_[jel] -> getRhsVector();              

                    //Disperse local contributions into the global matrix
                    for (int i=0; i<6; i++){
                        for (int j=0; j<6; j++){
                            if (fabs(Ajac(2*i  ,2*j  )) >= 1.e-15){
                                int dof_i = 3 * numNodesCoarse + 2 * connec(i);
                                int dof_j = 3 * numNodesCoarse + 2 * connec(j);
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    
                                                    &Ajac(2*i  ,2*j  ),
                                                    ADD_VALUES);
                            };
                            
                            if (fabs(Ajac(2*i+1,2*j  )) >= 1.e-15){
                                int dof_i = 3 * numNodesCoarse +
                                    2 * connec(i) + 1;
                                int dof_j = 3 * numNodesCoarse + 2 * connec(j);
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    
                                                    &Ajac(2*i+1,2*j  ),
                                                    ADD_VALUES);
                            };

                            if (fabs(Ajac(2*i  ,2*j+1)) >= 1.e-15){
                                int dof_i = 3 * numNodesCoarse + 2 * connec(i);
                                int dof_j = 3 * numNodesCoarse + 
                                    2 * connec(j) + 1;
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    
                                                    &Ajac(2*i  ,2*j+1),
                                                    ADD_VALUES);
                            };

                            if (fabs(Ajac(2*i+1,2*j+1)) >= 1.e-15){
                                int dof_i = 3 * numNodesCoarse + 
                                    2 * connec(i) + 1;
                                int dof_j = 3 * numNodesCoarse + 
                                    2 * connec(j) + 1;
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    
                                                    &Ajac(2*i+1,2*j+1),
                                                    ADD_VALUES);
                            };

                            //Matrix Q and Qt
                            if (fabs(Ajac(2*i  ,12+j)) >= 1.e-15){
                                int dof_i = 3 * numNodesCoarse + 2 * connec(i);
                                int dof_j = 3 * numNodesCoarse + 
                                    2 * numNodesFine + connec(j);
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,
                                                    &Ajac(2*i  ,12+j),
                                                    ADD_VALUES);
                            };

                            if (fabs(Ajac(12+j,2*i  )) >= 1.e-15){
                                int dof_i = 3 * numNodesCoarse + 2 * connec(i);
                                int dof_j = 3 * numNodesCoarse + 
                                    2 * numNodesFine + connec(j);
                                ierr = MatSetValues(A,1,&dof_j,1,&dof_i,  
                                                    &Ajac(12+j,2*i  ),
                                                    ADD_VALUES);
                            };

                            if (fabs(Ajac(2*i+1,12+j)) >= 1.e-15){
                                int dof_i = 3 * numNodesCoarse + 
                                    2 * connec(i) + 1;
                                int dof_j = 3 * numNodesCoarse + 
                                    2 * numNodesFine + connec(j);
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,
                                                    &Ajac(2*i+1,12+j),
                                                    ADD_VALUES);
                            };
                            
                            if (fabs(Ajac(12+j,2*i+1)) >= 1.e-15){
                                int dof_i = 3 * numNodesCoarse + 
                                    2 * connec(i) + 1;
                                int dof_j = 3 * numNodesCoarse + 
                                    2 * numNodesFine + connec(j);
                                ierr = MatSetValues(A,1,&dof_j,1,&dof_i,  
                                                    &Ajac(12+j,2*i+1),
                                                    ADD_VALUES);
                            };
                            
                            if (fabs(Ajac(12+i,12+j)) >= 1.e-15){
                                int dof_i = 3 * numNodesCoarse + 
                                    2 * numNodesFine + connec(i);
                                int dof_j = 3 * numNodesCoarse + 
                                    2 * numNodesFine + connec(j);
                                ierr = MatSetValues(A,1,&dof_j,1,&dof_i,  
                                                    &Ajac(12+i,12+j),
                                                    ADD_VALUES);
                            };
                        };     
                        ///Rhs vector
                        if (fabs(Rhs(2*i  )) >= 1.e-15){
                            int dof_i = 3 * numNodesCoarse + 2 * connec(i);
                            ierr = VecSetValues(b,1,&dof_i,&Rhs(2*i  ),
                                                ADD_VALUES);
                        };

                        if (fabs(Rhs(2*i+1)) >= 1.e-15){
                            int dof_i = 3 * numNodesCoarse + 2 * connec(i) + 1;
                            ierr = VecSetValues(b,1,&dof_i,&Rhs(2*i+1),
                                                ADD_VALUES);
                        };

                        if (fabs(Rhs(12+i)) >= 1.e-15){
                            int dof_i = 3 * numNodesCoarse + 2 * numNodesFine +
                                connec(i);
                            ierr = VecSetValues(b,1,&dof_i,&Rhs(12+i),
                                                ADD_VALUES);
                        };
                    }; 
                };                
            };
              

            //Lagrange Multipliers
            for (int l=0; l< numElemGlueZoneFine; l++){
                
                int jel = elementsGlueZoneFine_[l];
                
                if (domDecompFine.first[jel] == rank) {
                    
                    typename Elements::LocalMatrix Ajac, AjacAnt, AStab;
                    typename Elements::LocalVector Rhs, RhsStab;
                    typename Elements::Connectivity connec, connecC, connecL;
                    
                    connec = elementsFine_[jel] -> getConnectivity();
                    connecL = glueZoneFine_[l] -> getConnectivity();

                    // FINE MESH
                    
                    //Computes element matrix
                    elementsFine_[jel] -> getLagrangeMultipliersSameMesh();
                    
                    //Gets element matrice and rhs vectors
                    Ajac = - elementsFine_[jel] -> getLagrMultMatrix();
                    Rhs = - elementsFine_[jel] -> getRhsVector();
                    // AjacAnt = - elementsFine_[jel] -> getJacNRMatrix();
                    Ajac = trans(Ajac);

                    // elementsFine_[jel] -> getLMStabilizationSameMesh();

                    // AStab = - elementsFine_[jel] -> getJacNRMatrix();
                    // RhsStab = - elementsFine_[jel] -> getRhsVector();

                    //AjacAnt = trans(AjacAnt);
                    
                    typename Elements::LocalVector rhsLagMult,U_,lagStab;
                    
                    U_.clear();
                    rhsLagMult.clear();
                    lagStab.clear();
                    for (int i = 0; i < 6; i++){
                        U_(2*i  ) = nodesFine_[connec(i)] -> 
                            getLagrangeMultiplier(0);
                        U_(2*i+1) = nodesFine_[connec(i)] -> 
                            getLagrangeMultiplier(1);
                    };
                    
                    noalias(rhsLagMult) = -prod(trans(Ajac),U_);
                    // noalias(lagStab) = -prod(trans(AjacAnt),U_);
                   
                    
                    //Disperse local contributions into the global matrix
                    for (int i=0; i<6; i++){
                        for (int j=0; j<6; j++){
                            
                            if (fabs(Ajac(2*i  ,2*j  )) >= 1.e-15){
                                int d_i = 3 * numNodesCoarse + 3 * numNodesFine
                                    + 2 * connecL(i);
                                int d_j = 3 * numNodesCoarse + 2 * connec(j);
                                ierr = MatSetValues(A,1,&d_i,1,&d_j,
                                                    &Ajac(2*i  ,2*j  ),
                                                    ADD_VALUES);
                                ierr = MatSetValues(A,1,&d_j,1,&d_i,
                                                    &Ajac(2*i  ,2*j  ),
                                                    ADD_VALUES);
                            };
                            // ierr = MatSetValues(A,1,&d_i,1,&d_j,
                            //           &AStab(2*i  ,2*j  ),ADD_VALUES);

                            if (fabs(Ajac(2*i+1,2*j  )) >= 1.e-15){
                                int d_i = 3 * numNodesCoarse + 3 * numNodesFine
                                    + 2 * connecL(i) + 1;
                                int d_j = 3 * numNodesCoarse + 2 * connec(j);
                                ierr = MatSetValues(A,1,&d_i,1,&d_j,
                                                    &Ajac(2*i+1,2*j  ),
                                                    ADD_VALUES);
                                ierr = MatSetValues(A,1,&d_j,1,&d_i,
                                                    &Ajac(2*i+1,2*j  ),
                                                    ADD_VALUES);
                            };
                            // ierr = MatSetValues(A,1,&d_i,1,&d_j,
                            //                     &AStab(2*i+1,2*j  ),ADD_VALUES);

                            if (fabs(Ajac(2*i+1,2*j+1)) >= 1.e-15){
                                int d_i = 3 * numNodesCoarse + 3 * numNodesFine
                                    + 2 * connecL(i) + 1;
                                int d_j = 3 * numNodesCoarse + 2*connec(j) + 1;
                                ierr = MatSetValues(A,1,&d_i,1,&d_j,
                                                    &Ajac(2*i+1,2*j+1),
                                                    ADD_VALUES);
                                ierr = MatSetValues(A,1,&d_j,1,&d_i,
                                                    &Ajac(2*i+1,2*j+1),
                                                    ADD_VALUES);
                            };
                            // ierr = MatSetValues(A,1,&d_i,1,&d_j,
                            //                     &AStab(2*i+1,2*j+1),ADD_VALUES);

                            if (fabs(Ajac(2*i  ,2*j+1)) >= 1.e-15){
                                int d_i = 3 * numNodesCoarse + 3 * numNodesFine
                                    + 2 * connecL(i);
                                int d_j = 3 * numNodesCoarse + 2*connec(j) + 1;
                                ierr = MatSetValues(A,1,&d_i,1,&d_j,
                                                    &Ajac(2*i  ,2*j+1),
                                                    ADD_VALUES);
                                ierr = MatSetValues(A,1,&d_j,1,&d_i,
                                                    &Ajac(2*i  ,2*j+1),
                                                    ADD_VALUES);
                            };
                            // ierr = MatSetValues(A,1,&d_i,1,&d_j,
                            //                     &AStab(2*i  ,2*j+1),ADD_VALUES);
                            // if (fabs(AjacAnt(2*i  ,2*j  )) >= 1.e-15){
                            //     int d_i = 3 * numNodesCoarse + 3 * numNodesFine
                            //         + 2 * connecL(i);
                            //     int d_j = 3 * numNodesCoarse + 3 * numNodesFine
                            //         + 2 * connecL(j);
                            //     ierr = MatSetValues(A,1,&d_i,1,&d_j,
                            //                         &AjacAnt(2*i  ,2*j  )
                            //                         ,ADD_VALUES);
                            // };

                            // if (fabs(AjacAnt(2*i+1,2*j+1)) >= 1.e-15){
                            //     int d_i = 3 * numNodesCoarse + 3 * numNodesFine
                            //         + 2 * connecL(i) + 1;
                            //     int d_j = 3 * numNodesCoarse + 3 * numNodesFine
                            //         + 2 * connecL(j) + 1;
                            //     ierr = MatSetValues(A,1,&d_i,1,&d_j,
                            //                         &AjacAnt(2*i+1,2*j+1)
                            //                         ,ADD_VALUES);
                            // };
                        };
                        //Rhs vector

                        int dof_i = 3 * numNodesCoarse + 3 * numNodesFine +
                            2 * connecL(i);
                        ierr = VecSetValues(b,1,&dof_i,&Rhs(2*i  ),
                                            ADD_VALUES);
                        // ierr = VecSetValues(b,1,&dof_i,&RhsStab(2*i  ),
                        //                     ADD_VALUES);
                        // ierr = VecSetValues(b,1,&dof_i,&lagStab(2*i  ),
                        //                     ADD_VALUES);

                        dof_i = 3 * numNodesCoarse + 3 * numNodesFine + 
                            2 * connecL(i) + 1;
                        ierr = VecSetValues(b,1,&dof_i,&Rhs(2*i+1),
                                            ADD_VALUES);
                        // ierr = VecSetValues(b,1,&dof_i,&RhsStab(2*i+1),
                        //                     ADD_VALUES);
                        // ierr = VecSetValues(b,1,&dof_i,&lagStab(2*i+1),
                        //                     ADD_VALUES);

                        dof_i = 3 * numNodesCoarse + 2 * connec(i);
                        ierr = VecSetValues(b,1,&dof_i,&rhsLagMult(2*i  )
                                            ,ADD_VALUES);

                        dof_i = 3 * numNodesCoarse + 2 * connec(i) + 1;
                        ierr = VecSetValues(b,1,&dof_i,&rhsLagMult(2*i+1)
                                            ,ADD_VALUES);
                    };      

                    //COAESE MESH
                    
                    //Counts number of coarse mesh intersecting the 
                    //fine element
                    int numberIntPoints = elementsFine_[jel] -> 
                        getNumberOfIntegrationPoints();
                    int aux;
                    
                    std::vector<int> ele, diffElem;
                    ele.clear();
                    diffElem.clear();

                    ele.reserve(3);
                    for (int i=0; i<numberIntPoints; i++){
                        aux = elementsFine_[jel] -> 
                            getIntegPointCorrespondenceElement(i);
                        ele.push_back(aux);
                        //std::cout << "Num elem inters " << jel << " " << aux << std::endl;

                    };
                    
                    int numElemIntersect = 1;
                    int flag = 0;
                    diffElem.push_back(ele[0]);
                    
                    for (int i = 1; i<numberIntPoints; i++){
                        flag = 0;
                        for (int j = 0; j<numElemIntersect; j++){
                            if (ele[i] == diffElem[j]) {
                                break;
                            }else{
                                flag++;
                            };
                            if(flag == numElemIntersect){
                                numElemIntersect++;
                                diffElem.push_back(ele[i]);
                            };
                        };
                    };
                    //std::cout << "JEL " << jel << " " << numElemIntersect << std::endl;
                    //Compute the Lagrange Multiplier element matrix
                    for (int ielem = 0; ielem < numElemIntersect; ielem++){
                        
                        int iElemCoarse = diffElem[ielem];
                       
                        elementsFine_[jel] -> 
                            getLagrangeMultipliersDifferentMesh(iElemCoarse);
                        
                        //Computes element matrix
                        Ajac = elementsFine_[jel] -> getLagrMultMatrix();
                        Ajac = trans(Ajac);
                        // AjacAnt = elementsFine_[jel] -> getJacNRMatrix();
                        // AjacAnt = trans(AjacAnt);

                        connecC = elementsCoarse_[iElemCoarse] -> 
                            getConnectivity();

                        
                        typename Elements::LocalVector rhsLagMult,U_,lagStab;
                        
                        U_.clear();
                        for (int i = 0; i < 6; i++){
                            U_(2*i  ) = nodesFine_[connec(i)] -> 
                                getLagrangeMultiplier(0);
                            U_(2*i+1) = nodesFine_[connec(i)] -> 
                                getLagrangeMultiplier(1);
                        };
                        
                        noalias(rhsLagMult) = -prod(trans(Ajac),U_);
                        // noalias(lagStab) = -prod(trans(AjacAnt),U_);
                        
                        std::pair<Elements::LocalVector, 
                                  Elements::LocalMatrix>  lagMult;
                        
                        lagMult.first.clear();
                        lagMult.second.clear();
                         
                        lagMult = elementsCoarse_[iElemCoarse] -> 
                            getRhsVectorAndBoundaryConditions(Ajac);
                        
                        //Disperse local contribution into the global matrix
                        for (int i=0; i<6; i++){
                            for (int j=0; j<6; j++){

                                if (fabs(Ajac(2*i  ,2*j  )) >= 1.e-15){
                                    int dof_i = 3 * numNodesCoarse 
                                        + 3 * numNodesFine + 2 * connecL(i);
                                    int dof_j = 2 * connecC(j);
                                    ierr = MatSetValues(A,1,&dof_i,1,&dof_j,
                                                    &lagMult.second(2*i  ,2*j  )
                                                        ,ADD_VALUES);
                                    ierr = MatSetValues(A,1,&dof_j,1,&dof_i,
                                                    &lagMult.second(2*i  ,2*j  )
                                                        ,ADD_VALUES);
                                };

                                if (fabs(Ajac(2*i+1,2*j  )) >= 1.e-15){
                                    int dof_i = 3 * numNodesCoarse 
                                        + 3 * numNodesFine + 2 * connecL(i) + 1;
                                    int dof_j = 2 * connecC(j);
                                    ierr = MatSetValues(A,1,&dof_i,1,&dof_j,
                                                    &lagMult.second(2*i+1,2*j  )
                                                        ,ADD_VALUES);
                                    ierr = MatSetValues(A,1,&dof_j,1,&dof_i,
                                                    &lagMult.second(2*i+1,2*j  )
                                                        ,ADD_VALUES);
                                };

                                if (fabs(Ajac(2*i+1,2*j+1)) >= 1.e-15){
                                    int dof_i = 3 * numNodesCoarse 
                                        + 3 * numNodesFine + 2 * connecL(i) + 1;
                                    int dof_j = 2 * connecC(j) + 1;
                                    ierr = MatSetValues(A,1,&dof_i,1,&dof_j,
                                                    &lagMult.second(2*i+1,2*j+1)
                                                        ,ADD_VALUES);
                                    ierr = MatSetValues(A,1,&dof_j,1,&dof_i,
                                                    &lagMult.second(2*i+1,2*j+1)
                                                        ,ADD_VALUES);
                                };

                                if (fabs(Ajac(2*i  ,2*j+1)) >= 1.e-15){
                                    int dof_i = 3 * numNodesCoarse 
                                        + 3 * numNodesFine + 2 * connecL(i);
                                    int dof_j = 2 * connecC(j) + 1;
                                    ierr = MatSetValues(A,1,&dof_i,1,&dof_j,
                                                    &lagMult.second(2*i  ,2*j+1)
                                                        ,ADD_VALUES);
                                    ierr = MatSetValues(A,1,&dof_j,1,&dof_i,
                                                    &lagMult.second(2*i  ,2*j+1)
                                                        ,ADD_VALUES);
                                };

                                // if (fabs(AjacAnt(2*i  ,2*j  )) >= 1.e-15){
                                //     int d_i = 3 * numNodesCoarse + 
                                //         3 * numNodesFine + 2 * connecL(i);
                                //     int d_j = 3 * numNodesCoarse + 
                                //         3 * numNodesFine + 2 * connecL(j);
                                //     ierr = MatSetValues(A,1,&d_i,1,&d_j,
                                //                         &AjacAnt(2*i  ,2*j  ),
                                //                         ADD_VALUES);
                                // };
                                // if (fabs(AjacAnt(2*i+1,2*j+1)) >= 1.e-15){
                                //     int d_i = 3 * numNodesCoarse + 
                                //         3 * numNodesFine + 2 * connecL(i) + 1;
                                //     int d_j = 3 * numNodesCoarse + 
                                //         3 * numNodesFine + 2 * connecL(j) + 1;
                                //     ierr = MatSetValues(A,1,&d_i,1,&d_j,
                                //                         &AjacAnt(2*i+1,2*j+1),
                                //                         ADD_VALUES);
                                // };
                            };
                            //Rhs vector

                            int d_i = 3 * numNodesCoarse 
                                + 3 * numNodesFine + 2 * connecL(i);
                            ierr = VecSetValues(b,1,&d_i,
                                                &lagMult.first(2*i  ),
                                                ADD_VALUES);
                            // ierr = VecSetValues(b,1,&d_i,
                            //                     &lagStab(2*i  ),
                            //                     ADD_VALUES);

                            d_i = 3 * numNodesCoarse 
                                + 3 * numNodesFine + 2 * connecL(i) + 1;
                            ierr = VecSetValues(b,1,&d_i,
                                                &lagMult.first(2*i+1),
                                                ADD_VALUES);
                            // ierr = VecSetValues(b,1,&d_i,
                            //                     &lagStab(2*i+1),
                            //                     ADD_VALUES);

                            int dof_i = 2 * connecC(i);
                            ierr = VecSetValues(b,1,&dof_i,
                                                &rhsLagMult(2*i  ),
                                                ADD_VALUES);

                            dof_i = 2 * connecC(i) + 1;
                            ierr = VecSetValues(b,1,&dof_i,
                                                &rhsLagMult(2*i+1),
                                                ADD_VALUES);
                        };                                 
                    }; //Number of intersections
                }; // if element belongs to the glue zone
            }; // Glue zone


            //------------------------------------------------------------------
            //-------------------END OF LINEAR SYSTEM ASSEMBLY------------------
            //------------------------------------------------------------------
            //std::cout << "Enter PETSc " << rank << std::endl;
            
            //Assemble matrices and vectors
            ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
            ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
            
            ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
            ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
            
            //ierr = MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
            // ierr = VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
            
            //Create KSP context to solve the linear system
            ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
            
            ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
            



            // ierr = MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE,0, NULL, &nullsp);
            // ierr = MatSetNullSpace(A, nullsp);
            // ierr = MatNullSpaceDestroy(&nullsp);




            // ierr = KSPSetTolerances(ksp,1.e-7,1.e-10,PETSC_DEFAULT,
            //                         1000);CHKERRQ(ierr);
            
            // //ierr = KSPGMRESSetRestart(ksp, 10); CHKERRQ(ierr);
            
            // ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
            
            // ierr = PCSetType(pc,PCNONE);CHKERRQ(ierr);
            
            // //ierr = KSPSetPCSide(ksp, PC_RIGHT);
            // ierr = KSPSetType(ksp,KSPLSQR); CHKERRQ(ierr);
            
            // ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
            // // ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);



            
#if defined(PETSC_HAVE_MUMPS)
            ierr = KSPSetType(ksp,KSPPREONLY);
            ierr = KSPGetPC(ksp,&pc);
            ierr = PCSetType(pc, PCLU);
#endif
            
            ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
            ierr = KSPSetUp(ksp);
            




           // ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);

            ierr = KSPSolve(ksp,b,u);CHKERRQ(ierr);
            
            ierr = KSPGetTotalIterations(ksp, &iterations);
            
            //if (rank == 0)std::cout << "GMRES Iterations = " << iterations << std::endl;
            
            //ierr = VecView(u,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
            
            //Gathers the solution vector to the master process
            ierr = VecScatterCreateToAll(u, &ctx, &All);CHKERRQ(ierr);
            
            ierr = VecScatterBegin(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);
            CHKERRQ(ierr);
            
            ierr = VecScatterEnd(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);
            CHKERRQ(ierr);
            
            ierr = VecScatterDestroy(&ctx);CHKERRQ(ierr);
            
            //Updates nodal values
            double u_[2];
            double norm = 0.;
            double p_;
            Ione = 1;
            
            for (int i = 0; i < numNodesCoarse; ++i){
                Ii = 2 * i;
                ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                u_[0] = val;
                norm += val*val;
                nodesCoarse_[i] -> incrementVelocity(0,u_[0]);
                
                Ii = 2 * i + 1;
                ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                u_[1] = val;
                norm += val*val;
                nodesCoarse_[i] -> incrementVelocity(1,u_[1]);
            };
            
            for (int i = 0; i<numNodesCoarse; i++){
                Ii = 2 * numNodesCoarse + i;
                ierr = VecGetValues(All,Ione,&Ii,&val);CHKERRQ(ierr);
                p_ = val;
                nodesCoarse_[i] -> incrementPressure(p_);
            };
            
            for (int i = 0; i < numNodesFine; ++i){
                Ii = 3 * numNodesCoarse + 2 * i;
                ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                u_[0] = val;
                norm += val*val;
                nodesFine_[i] -> incrementVelocity(0,u_[0]);
                
                Ii = 3 * numNodesCoarse + 2 * i + 1;
                ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                u_[1] = val;
                norm += val*val;
                nodesFine_[i] -> incrementVelocity(1,u_[1]);
            };
            
            for (int i = 0; i<numNodesFine; i++){
                Ii = 3 * numNodesCoarse + 2 * numNodesFine + i;
                ierr = VecGetValues(All,Ione,&Ii,&val);CHKERRQ(ierr);
                p_ = val;
                nodesFine_[i] -> incrementPressure(p_);
            };
            
            for (int i = 0; i < numNodesGlueZoneFine; ++i){
                Ii = 3 * numNodesCoarse + 3 * numNodesFine + 2 * i;
                ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                u_[0] = val;
                nodesFine_[nodesGlueZoneFine_[i]] -> incrementLagrangeMultiplier(0,u_[0]);
                
                Ii = 3 * numNodesCoarse + 3 * numNodesFine + 2 * i + 1;
                ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                u_[1] = val;
                nodesFine_[nodesGlueZoneFine_[i]] -> incrementLagrangeMultiplier(1,u_[1]);
                
                // std::cout << "LAG M " << u_[0] << " " << u_[1] << std::endl;
            };
            
            //Computes the solution vector norm
            ierr = VecNorm(u,NORM_2,&val);CHKERRQ(ierr);
            
            boost::posix_time::ptime t2 =
                boost::posix_time::microsec_clock::local_time();
            
            if(rank == 0){
                boost::posix_time::time_duration diff = t2 - t1;
                
                std::cout<<"Iteration = " << inewton << " (" << iterations <<  
                    ")  Du Norm = " << std::scientific << sqrt(norm) 
                         << " " << val << 
                    "  Time (s) = " << std::fixed << 
                    diff.total_milliseconds()/1000. << std::endl;
            };
            
            ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
            ierr = VecDestroy(&b); CHKERRQ(ierr);
            ierr = VecDestroy(&u); CHKERRQ(ierr);
            ierr = VecDestroy(&All); CHKERRQ(ierr);
            ierr = MatDestroy(&A); CHKERRQ(ierr);
            
            if(val <= tolerance){
                break;            
            };          

            //Updates SUPG Parameter
            // for (int i = 0; i < numElemFine; i++){
            //     elementsFine_[i] -> getParameterSUPG();
            // };   
            // for (int i = 0; i < numElemCoarse; i++){
            //     elementsCoarse_[i] -> getParameterSUPG();
            // };

        };

        //Compute real velocity
        
        QuadShapeFunction<2>                       shapeQuad;
        typename QuadShapeFunction<2>::Values      phi_;
        
        for (int i = 0; i<numNodesFine; i++){
            nodesFine_[i] -> setVelocityArlequin(0,nodesFine_[i] -> 
                                                 getVelocity(0));
            nodesFine_[i] -> setVelocityArlequin(1,nodesFine_[i] -> 
                                                 getVelocity(1));
        };
        
        for (int i = 0; i<numNodesGlueZoneFine; i++){
            typename Nodes::VecLocD xsi;
            typename Quadrature::NodalValuesQuad u_coarse, v_coarse;
            typename Elements::Connectivity connecCoarse;
            
            double u = 0.;
            double v = 0.;
            
            int elCoarse = nodesFine_[nodesGlueZoneFine_[i]] -> 
                getNodalElemCorrespondence();
            xsi = nodesFine_[nodesGlueZoneFine_[i]] -> 
                getNodalXsiCorrespondence();
            
            connecCoarse = elementsCoarse_[elCoarse] -> getConnectivity();
            
            for (int j=0; j<6; j++){
                u_coarse(j) = nodesCoarse_[connecCoarse(j)] -> getVelocity(0);
                v_coarse(j) = nodesCoarse_[connecCoarse(j)] -> getVelocity(1);
            };
            
            shapeQuad.evaluate(xsi,phi_);
            
            for (int j=0; j<6; j++){
                u += u_coarse(j) * phi_(j);
                v += v_coarse(j) * phi_(j);
            };
            
            double wFunc = nodesFine_[nodesGlueZoneFine_[i]] -> 
                getWeightFunction();
            
            double u_int = nodesFine_[nodesGlueZoneFine_[i]] -> 
                getVelocity(0) * wFunc + u * (1. - wFunc);
            double v_int = nodesFine_[nodesGlueZoneFine_[i]] -> 
                getVelocity(1) * wFunc + v * (1. - wFunc);
            
            nodesFine_[nodesGlueZoneFine_[i]] -> setVelocityArlequin(0,u_int);
            nodesFine_[nodesGlueZoneFine_[i]] -> setVelocityArlequin(1,v_int);
            
        };
        
        for (int i=0; i<numNodesCoarse; i++){
            nodesCoarse_[i] -> 
                setVelocityArlequin(0,nodesCoarse_[i] -> getVelocity(0));
            nodesCoarse_[i] -> 
                setVelocityArlequin(1,nodesCoarse_[i] -> getVelocity(1));
        };
        
        //Compute real pressure
        
        // LinShapeFunction<2>                       shapeLin;
        // typename LinShapeFunction<2>::Values      phip_;
        
        // for (int i=0; i<numNodesPFine; i++){
        //     nodesPFine_[i] -> 
        //         setPressureArlequin(nodesPFine_[i] -> getPressure());
        // };
        
        // for (int jel = 0; jel<numElemGlueZoneFine; jel++){
            
        //     typename Nodes::VecLocD xsi;
        //     typename Quadrature::NodalValuesLin p_coarse;
        //     typename Elements::ConnectivityP connecPCoarse, connecPFine;
        //     typename Elements::Connectivity connecCoarse, connec;
            
        //     connec = elementsFine_[elementsGlueZoneFine_[jel]] -> 
        //         getConnectivity();
        //     connecPFine = elementsFine_[elementsGlueZoneFine_[jel]] -> 
        //         getConnectivityP();
            
        //     for (int i=0; i<3; i++){
                
        //         double p = 0.;
                
        //         int elCoarse = nodesFine_[connec(i)] ->
        //             getNodalElemCorrespondence();
        //         xsi = nodesFine_[connec(i)] -> getNodalXsiCorrespondence();
                
        //         connecCoarse = elementsCoarse_[elCoarse] -> getConnectivity();
        //         connecPCoarse = elementsCoarse_[elCoarse] -> getConnectivityP();
                
        //         for (int j=0; j<3; j++){
        //             p_coarse(j) = nodesPCoarse_[connecPCoarse(j)] -> 
        //                 getPressure();
        //         };
                
        //         shapeLin.evaluate(xsi,phip_);
                
        //         for (int j=0; j<3; j++){
        //             p += p_coarse(j) * phip_(j);
        //         };
                
        //         //double wFunc = nodesFine_[connec(i)] -> getWeightFunction();
                
        //         typename Nodes::VecLocD x;

        //         double p_int = 0.;
       
        //         if (fabs(nodesPFine_[connecPFine(i)] -> getPressure()) > 0){
        //             p_int = nodesPFine_[connecPFine(i)] -> getPressure();
        //         }else{
        //             p_int = nodesPFine_[connecPFine(i)] -> getPressure() + p;
        //         };

        //         nodesPFine_[connecPFine(i)] -> setPressureArlequin(p_int);
        //     };
        // };
        
        // for (int i=0; i<numNodesPCoarse; i++){
        //     nodesPCoarse_[i] -> 
        //         setPressureArlequin(nodesPCoarse_[i] -> getPressure());
        // };
        
        
        if (rank == 0) {
            //Printing results
            printVelocity(iTimeStep);
        };

    };
        
    return 0;

};

//------------------------------------------------------------------------------
//----------------COMPUTE ARLEQUIN COUPLED NAVIER-STOKES PROBLEM----------------
//------------------------------------------------------------------------------
template<>
int Arlequin<2>::solveArlequinProblemCoarse(int iterNumber, double tolerance,
                                            int problem_type, 
                                            int time_dependency){

    // Mat               A;
    // Vec               b, u, All;
    // PetscErrorCode    ierr;
    // PetscInt          Istart, Iend, Ii, Ione, iterations;
    // KSP               ksp; 
    // PC                pc;
    // VecScatter        ctx;
    // PetscScalar       val;
    // //MatNullSpace      nullsp;


    // if ((problem_type < 1) || (problem_type > 2)){
    //     std::cout << "WRONG PROBLEM TYPE." << std::endl;
    //     return 0;
    // };

    // if (time_dependency == 0) numTimeSteps = 1;

    
    // //Construct the glue zone based on some defined criterion
    // setCouplingZoneCoarse();

    // //Computes the Weight function for all the finite elements
    // setWeightFunction(4.);
   
    // //Computes the Nodal correspondence between fine nodes and coarse elements
    // setNodalCorrespondenceFine();

    // //Construct the free zone in coarse model
    // setFreeZone();

    // // for (int jel = 0; jel < numElemGlueZoneFine; jel++){
    // //     int iel = elementsGlueZoneFine_[jel];
    // //     elementsFine_[iel] -> setCompressibility();
    // // };

    // MPI_Comm_rank(PETSC_COMM_WORLD, &rank);      

    // // Computes the system size
    // int sysSize = 3 * numNodesCoarse + 
    //     + 3 * numNodesFine + 2 * numNodesCoarse;
    
    // for (int iTimeStep = 0; iTimeStep < numTimeSteps; iTimeStep++){
        
    //     if (rank == 0) {std::cout << "------------------------- TIME STEP = "
    //                               << iTimeStep << " -------------------------"
    //                               << std::endl;}
        
    //     //Updates velocity and acceleration
    //     for (int i = 0; i < numNodesCoarse; i++){
    //         double accel[2], u[2], uprev[2];
            
    //         //Compute acceleration
    //         u[0] = nodesCoarse_[i] -> getVelocity(0);
    //         u[1] = nodesCoarse_[i] -> getVelocity(1);
            
    //         uprev[0] = nodesCoarse_[i] -> getPreviousVelocity(0);
    //         uprev[1] = nodesCoarse_[i] -> getPreviousVelocity(1);
            
    //         accel[0] = (u[0] - uprev[0]) / dTime;
    //         accel[1] = (u[1] - uprev[1]) / dTime;
            
    //         nodesCoarse_[i] -> setAcceleration(accel);
            
    //         //Updates velocity
    //         nodesCoarse_[i] -> setPreviousVelocity(u);
            
    //         //Computes predictor
    //         u[0] += dTime * accel[0];
    //         u[1] += dTime * accel[1];
            
    //         nodesCoarse_[i] -> setVelocity(u);
            
    //         //Updates acceleration
    //         nodesCoarse_[i] -> setPreviousAcceleration(accel);
            
    //     };
    //     for (int i = 0; i < numNodesFine; i++){
    //         double accel[2], u[2], uprev[2];
            
    //         //Compute acceleration
    //         u[0] = nodesFine_[i] -> getVelocity(0);
    //         u[1] = nodesFine_[i] -> getVelocity(1);
            
    //         uprev[0] = nodesFine_[i] -> getPreviousVelocity(0);
    //         uprev[1] = nodesFine_[i] -> getPreviousVelocity(1);
            
    //         accel[0] = (u[0] - uprev[0]) / dTime;
    //         accel[1] = (u[1] - uprev[1]) / dTime;
            
    //         nodesFine_[i] -> setAcceleration(accel);
            
    //         //Updates velocity
    //         nodesFine_[i] -> setPreviousVelocity(u);
            
    //         //Computes predictor
    //         u[0] += dTime * accel[0];
    //         u[1] += dTime * accel[1];
            
    //         nodesFine_[i] -> setVelocity(u);
            
    //         //Updates acceleration
    //         nodesFine_[i] -> setPreviousAcceleration(accel);
            
    //     };

        //STARTS NEWTON-RAPHSON
//         for (int inewton = 0; inewton < iterNumber; inewton++){
            
//             // if (inewton < 10){
//             //     for (int i = 0; i < numElemFine; i++){
//             //         double visc = 0.001*(iTimeStep+1)/10.;
//             //         elementsFine_[i] -> setViscosity(visc);
//             //     };
//             //     for (int i = 0; i < numElemCoarse; i++){
//             //         double visc = 0.001*(iTimeStep+1)/10.;
//             //         elementsCoarse_[i] -> setViscosity(visc);
//             //     };
//             // }else{
//             //     for (int i = 0; i < numElemFine; i++){
//             //         double visc = 0.001;
//             //         elementsFine_[i] -> setViscosity(visc);
//             //     };
//             //     for (int i = 0; i < numElemCoarse; i++){
//             //         double visc = 0.001;
//             //         elementsCoarse_[i] -> setViscosity(visc);
//             //     };
//             // };


//             boost::posix_time::ptime t1 =
//                 boost::posix_time::microsec_clock::local_time();
            
//             // Preallocates the matrix
//             ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
//                                 sysSize, sysSize, 160, NULL, 300, NULL, &A); 
//             CHKERRQ(ierr);
            
//             // Divides the matrix between the processes
//             ierr = MatGetOwnershipRange(A, &Istart, &Iend);CHKERRQ(ierr);
            
//             //Create PETSc vectors
//             ierr = VecCreate(PETSC_COMM_WORLD, &b); CHKERRQ(ierr);
//             ierr = VecSetSizes(b, PETSC_DECIDE, sysSize); CHKERRQ(ierr);
//             ierr = VecSetFromOptions(b); CHKERRQ(ierr); 
//             ierr = VecDuplicate(b, &u); CHKERRQ(ierr);
//             ierr = VecDuplicate(b, &All); CHKERRQ(ierr);
                        
//             //------------------------------------------------------------------
//             //------------------BEGIN OF LINEAR SYSTEM ASSEMBLY-----------------
//             //------------------------------------------------------------------


//             // for (int i=0; i<sysSize; i++){
//             //     int dof_i = i;
//             //     double one = 1.e-15;
//             //     ierr = MatSetValues(A,1,&dof_i,1,&dof_i,&one,ADD_VALUES);
//             // };


//             //Coarse mesh
            
//             for (int jel = 0; jel < numElemCoarse; jel++){   
                
//                 if (domDecompCoarse.first[jel] == rank) {
                    
//                     //Compute Element matrix
//                     if (time_dependency == 0){ 
//                         if (problem_type == 1){
//                             elementsCoarse_[jel] -> getSteadyStokes();
//                         }else{
//                             if (inewton == 0){
//                                 elementsCoarse_[jel] -> getSteadyStokes();
//                             }else{
//                                 elementsCoarse_[jel] -> getSteadyNavierStokes();
//                             };
//                         };
//                     }else{
//                         if (problem_type == 1){
//                             elementsCoarse_[jel] -> getTransientStokes();
//                         }else{
//                             elementsCoarse_[jel] -> getTransientNavierStokes();
//                         };
//                     };

//                     typename Elements::LocalMatrix Ajac;
//                     typename Elements::LocalVector Rhs;
//                     typename Elements::Connectivity connec;
//                     typename Elements::ConnectivityP connecP;
                    
//                     //Gets element connectivity, jacobian and rhs 
//                     connec = elementsCoarse_[jel] -> getConnectivity();
//                     connecP = elementsCoarse_[jel] -> getConnectivityP();
//                     Ajac = elementsCoarse_[jel] -> getJacNRMatrix();
//                     Rhs = elementsCoarse_[jel] -> getRhsVector();
                                      
//                     //Disperse local contributions into the global matrix
//                     for (int i=0; i<6; i++){
//                         for (int j=0; j<6; j++){
//                             if (fabs(Ajac(2*i  ,2*j  )) >= 1.e-15){
//                                 int dof_i = 2 * connec(i);
//                                 int dof_j = 2 * connec(j);
//                                 ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    
//                                                 &Ajac(2*i  ,2*j  ),ADD_VALUES);
//                             };
//                             if (fabs(Ajac(2*i+1,2*j  )) >= 1.e-15){
//                                 int dof_i = 2 * connec(i) + 1;
//                                 int dof_j = 2 * connec(j);
//                                 ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    
//                                                 &Ajac(2*i+1,2*j  ),ADD_VALUES);
//                             };
//                             if (fabs(Ajac(2*i  ,2*j+1)) >= 1.e-15){
//                                 int dof_i = 2 * connec(i);
//                                 int dof_j = 2 * connec(j) + 1;
//                                 ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    
//                                                 &Ajac(2*i  ,2*j+1),ADD_VALUES);
//                             };
//                             if (fabs(Ajac(2*i+1,2*j+1)) >= 1.e-15){
//                                 int dof_i = 2 * connec(i) + 1;
//                                 int dof_j = 2 * connec(j) + 1;
//                                 ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    
//                                                 &Ajac(2*i+1,2*j+1),ADD_VALUES);
//                             };
//                         };    
//                         //Matrix Q and Qt
//                         for (int j=0; j<3; j++){
//                             if (fabs(Ajac(2*i  ,12+j)) >= 1.e-15){
//                                 int dof_i = 2 * connec(i);
//                                 int dof_j = 2 * numNodesCoarse + connecP(j);
//                                 ierr = MatSetValues(A,1,&dof_i,1,&dof_j,
//                                                 &Ajac(2*i  ,12+j),ADD_VALUES);
//                                 ierr = MatSetValues(A,1,&dof_j,1,&dof_i,    
//                                                 &Ajac(12+j,2*i  ),ADD_VALUES);
//                             };
//                             if (fabs(Ajac(2*i+1,12+j)) >= 1.e-15){
//                                 int dof_i = 2 * connec(i) + 1;
//                                 int dof_j = 2 * numNodesCoarse + connecP(j);
//                                 ierr = MatSetValues(A,1,&dof_i,1,&dof_j,
//                                                 &Ajac(2*i+1,12+j),ADD_VALUES);
//                                 ierr = MatSetValues(A,1,&dof_j,1,&dof_i,    
//                                                 &Ajac(12+j,2*i+1),ADD_VALUES);
//                             };
//                         };            
       
//                         //Rhs vector
//                         if (fabs(Rhs(2*i  )) >= 1.e-15){
//                             int dof_i = 2 * connec(i);
//                             ierr = VecSetValues(b,1,&dof_i,&Rhs(2*i  ),
//                                                 ADD_VALUES);
//                         };
//                         if (fabs(Rhs(2*i+1)) >= 1.e-15){
//                             int dof_i = 2 * connec(i) + 1;
//                             ierr = VecSetValues(b,1,&dof_i,&Rhs(2*i+1),
//                                                 ADD_VALUES);
//                         };
//                     };
//                     for (int i=0; i<3; i++){
//                         if (fabs(Rhs(12+i)) >= 1.e-15){
//                             int dof_i = 2 * numNodesCoarse + connecP(i);
//                             ierr = VecSetValues(b,1,&dof_i,&Rhs(12+i),
//                                                 ADD_VALUES);
//                         };
//                     };
//                 };
//             };
            
//             //Fine mesh
//             for (int jel = 0; jel < numElemFine; jel++){   
                
//                 if (domDecompFine.first[jel] == rank) {
                    
//                     //Compute Element matrix
//                     if (time_dependency == 0) {
//                         if (problem_type == 1){
//                             elementsFine_[jel] -> getSteadyStokes();
//                         }else{
//                             if (inewton == 0){
//                                 elementsFine_[jel] -> getSteadyStokes();
//                             }else{
//                                 elementsFine_[jel] -> getSteadyNavierStokes();
//                             };
//                         };
//                     }else{
//                         if (problem_type == 1){
//                             elementsFine_[jel] -> getTransientStokes();
//                         }else{
//                             elementsFine_[jel] -> getTransientNavierStokes();
//                         };
//                     };

//                     typename Elements::LocalMatrix Ajac;
//                     typename Elements::LocalVector Rhs;
//                     typename Elements::Connectivity connec;
//                     typename Elements::ConnectivityP connecP;
                    
//                     //Gets element connectivity, jacobian and rhs 
//                     connec = elementsFine_[jel] -> getConnectivity();
//                     connecP = elementsFine_[jel] -> getConnectivityP();
//                     Ajac = elementsFine_[jel] -> getJacNRMatrix();
//                     Rhs = elementsFine_[jel] -> getRhsVector();
                    
//                     //Disperse local contributions into the global matrix
//                     for (int i=0; i<6; i++){
//                         for (int j=0; j<6; j++){
//                             if (fabs(Ajac(2*i  ,2*j  )) >= 1.e-15){
//                                 int dof_i = 2 * numNodesCoarse + 
//                                     numNodesPCoarse + 2 * connec(i);
//                                 int dof_j = 2 * numNodesCoarse + 
//                                     numNodesPCoarse + 2 * connec(j);
//                                 ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    
//                                                 &Ajac(2*i  ,2*j  ),ADD_VALUES);
//                             };
//                             if (fabs(Ajac(2*i+1,2*j  )) >= 1.e-15){
//                                 int dof_i = 2 * numNodesCoarse + 
//                                     numNodesPCoarse + 2 * connec(i) + 1;
//                                 int dof_j = 2 * numNodesCoarse + 
//                                     numNodesPCoarse + 2 * connec(j);
//                                 ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    
//                                                 &Ajac(2*i+1,2*j  ),ADD_VALUES);
//                             };
//                             if (fabs(Ajac(2*i  ,2*j+1)) >= 1.e-15){
//                                 int dof_i = 2 * numNodesCoarse + 
//                                     numNodesPCoarse + 2 * connec(i);
//                                 int dof_j = 2 * numNodesCoarse + 
//                                     numNodesPCoarse + 2 * connec(j) + 1;
//                                 ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    
//                                                 &Ajac(2*i  ,2*j+1),ADD_VALUES);
//                             };
//                             if (fabs(Ajac(2*i+1,2*j+1)) >= 1.e-15){
//                                 int dof_i = 2 * numNodesCoarse + 
//                                     numNodesPCoarse + 2 * connec(i) + 1;
//                                 int dof_j = 2 * numNodesCoarse + 
//                                     numNodesPCoarse + 2 * connec(j) + 1;
//                                 ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    
//                                                 &Ajac(2*i+1,2*j+1),ADD_VALUES);
//                             };
//                         };
//                         //Matrix Q and Qt
//                         for (int j=0; j<3; j++){
//                             if (fabs(Ajac(2*i  ,12+j)) >= 1.e-15){
//                                 int dof_i = 2 * numNodesCoarse + 
//                                     numNodesPCoarse + 2 * connec(i);
//                                 int dof_j = 2 * numNodesCoarse +
//                                     numNodesPCoarse + 2 * numNodesFine + 
//                                     connecP(j);
//                                 ierr = MatSetValues(A,1,&dof_i,1,&dof_j,
//                                                 &Ajac(2*i  ,12+j),ADD_VALUES);
//                                 ierr = MatSetValues(A,1,&dof_j,1,&dof_i,  
//                                                 &Ajac(12+j,2*i  ),ADD_VALUES);
//                             };
//                             if (fabs(Ajac(2*i+1,12+j)) >= 1.e-15){
//                                 int dof_i = 2 * numNodesCoarse + 
//                                     numNodesPCoarse + 2 * connec(i)+1;
//                                 int dof_j = 2 * numNodesCoarse + 
//                                     numNodesPCoarse + 2 * numNodesFine + 
//                                     connecP(j);
//                                 ierr = MatSetValues(A,1,&dof_i,1,&dof_j,
//                                                 &Ajac(2*i+1,12+j),ADD_VALUES);
//                                 ierr = MatSetValues(A,1,&dof_j,1,&dof_i,  
//                                                 &Ajac(12+j,2*i+1),ADD_VALUES);
//                             };
//                         };          
//                         //Rhs vector
//                         if (fabs(Rhs(2*i  )) >= 1.e-15){
//                             int dof_i = 2 * numNodesCoarse + numNodesPCoarse + 
//                                 2 * connec(i);
//                             ierr = VecSetValues(b,1,&dof_i,&Rhs(2*i  ),
//                                                 ADD_VALUES);
//                         };
//                         if (fabs(Rhs(2*i+1)) >= 1.e-15){
//                             int dof_i = 2 * numNodesCoarse + numNodesPCoarse + 
//                                 2 * connec(i) + 1;
//                             ierr = VecSetValues(b,1,&dof_i,&Rhs(2*i+1),
//                                                 ADD_VALUES);
//                         };
//                     };
//                     for (int i=0; i<3; i++){
//                         if (fabs(Rhs(12+i)) >= 1.e-15){
//                             int dof_i = 2 * numNodesCoarse + numNodesPCoarse + 
//                                 2 * numNodesFine + connecP(i);
//                             ierr = VecSetValues(b,1,&dof_i,&Rhs(12+i),
//                                                 ADD_VALUES);
//                         };
//                     };                                        
                    
//                 };
//             };
               
//             //Lagrange Multipliers
//             for (int l=0; l< numElemGlueZoneCoarse; l++){
                
//                 int jel = elementsGlueZoneCoarse_[l];

//                 if (domDecompCoarse.first[jel] == rank) {
                    
//                     typename Elements::LocalMatrix Ajac, AjacAnt;
//                     typename Elements::LocalVector Rhs;
//                     typename Elements::Connectivity connec, connecC;
                    
//                     connecC = elementsCoarse_[jel] -> getConnectivity();
                    
//                     // COARSE MESH
                    
//                     //Computes element matrix
//                     elementsCoarse_[jel] -> getLagrangeMultipliersSameMesh();
                    
//                     //Gets element matrice and rhs vectors
//                     Ajac = - elementsCoarse_[jel] -> getJacNRMatrix();
//                     Rhs = - elementsCoarse_[jel] -> getRhsVector();
//                     AjacAnt = - elementsCoarse_[jel] -> getLagrMultMatrix();

//                     typename Elements::LocalVector rhsLagMult,U_;
                    
//                     U_.clear();
//                     for (int i = 0; i < 6; i++){
//                         U_(2*i  ) = nodesCoarse_[connecC(i)] -> 
//                             getLagrangeMultiplier(0);
//                         U_(2*i+1) = nodesCoarse_[connecC(i)] -> 
//                             getLagrangeMultiplier(1);
//                     };
                    
//                     noalias(rhsLagMult) = - prod(trans(Ajac),U_);
                    
//                     //Disperse local contributions into the global matrix
//                     for (int i=0; i<6; i++){
//                         for (int j=0; j<6; j++){
//                             if (fabs(Ajac(2*i  ,2*j  )) >= 1.e-15){
//                                 int d_i = 2 * numNodesCoarse + numNodesPCoarse +
//                                     2 * numNodesFine + numNodesPFine + 
//                                     2 * connecC(i);
//                                 int d_j = 2 * connecC(j);
//                                 ierr = MatSetValues(A,1,&d_i,1,&d_j,
//                                                 &Ajac(2*i  ,2*j  ),ADD_VALUES);
//                                 ierr = MatSetValues(A,1,&d_j,1,&d_i,
//                                                 &Ajac(2*i  ,2*j  ),ADD_VALUES);
//                             };
//                             if (fabs(Ajac(2*i+1,2*j  )) >= 1.e-15){
//                                 int d_i = 2 * numNodesCoarse + numNodesPCoarse +
//                                     2 * numNodesFine + numNodesPFine + 
//                                     2 * connecC(i) + 1;
//                                 int d_j = 2 * connecC(j);
//                                 ierr = MatSetValues(A,1,&d_i,1,&d_j,
//                                                 &Ajac(2*i+1,2*j  ),ADD_VALUES);
//                                 ierr = MatSetValues(A,1,&d_j,1,&d_i,
//                                                 &Ajac(2*i+1,2*j  ),ADD_VALUES);
//                             };
//                             if (fabs(Ajac(2*i+1,2*j+1)) >= 1.e-15){
//                                 int d_i = 2 * numNodesCoarse + numNodesPCoarse +
//                                     2 * numNodesFine + numNodesPFine + 
//                                     2 * connecC(i) + 1;
//                                 int d_j = 2 * connecC(j) + 1;  
//                                 ierr = MatSetValues(A,1,&d_i,1,&d_j,
//                                                 &Ajac(2*i+1,2*j+1),ADD_VALUES);
//                                 ierr = MatSetValues(A,1,&d_j,1,&d_i,
//                                                 &Ajac(2*i+1,2*j+1),ADD_VALUES);
//                             };
//                             if (fabs(Ajac(2*i  ,2*j+1)) >= 1.e-15){
//                                 int d_i = 2 * numNodesCoarse + numNodesPCoarse +
//                                     2 * numNodesFine + numNodesPFine + 
//                                     2 * connecC(i);
//                                 int d_j = 2 * connecC(j) + 1;
//                                 ierr = MatSetValues(A,1,&d_i,1,&d_j,
//                                                 &Ajac(2*i  ,2*j+1),ADD_VALUES);
//                                 ierr = MatSetValues(A,1,&d_j,1,&d_i,
//                                                 &Ajac(2*i  ,2*j+1),ADD_VALUES);
//                             };
                            
//                         };
//                         // Rhs vector
//                         if (fabs(Rhs(2*i  )) >= 1.e-15){
//                             int dof_i = 2 * numNodesCoarse + numNodesPCoarse + 
//                                 2 * numNodesFine + numNodesPFine + 
//                                 2 * connecC(i);
//                             ierr = VecSetValues(b,1,&dof_i,&Rhs(2*i  ),
//                                                 ADD_VALUES);
//                         };
//                         if (fabs(Rhs(2*i+1)) >= 1.e-15){
//                             int dof_i = 2 * numNodesCoarse + numNodesPCoarse + 
//                                 2 * numNodesFine + numNodesPFine + 
//                                 2 * connecC(i) + 1;
//                             ierr = VecSetValues(b,1,&dof_i,&Rhs(2*i+1),
//                                                 ADD_VALUES);
//                         };
                        
//                         if (fabs(rhsLagMult(2*i  )) >= 1.e-15){
//                             int dof_i = 2 * connecC(i);
//                             ierr = VecSetValues(b,1,&dof_i,&rhsLagMult(2*i  )
//                                                 ,ADD_VALUES);
//                         };
//                         if (fabs(rhsLagMult(2*i+1)) >= 1.e-15){
//                             int dof_i = 2 * connecC(i) + 1;
//                             ierr = VecSetValues(b,1,&dof_i,&rhsLagMult(2*i+1)
//                                                 ,ADD_VALUES);
//                         };
//                     };      
//                 };
//             };
                    
//             //Lagrange Multipliers
//             for (int l=0; l< numElemGlueZoneFine; l++){
                
//                 int jel = elementsGlueZoneFine_[l];
                
//                 if (domDecompFine.first[jel] == rank) {

//                     //FINE MESH

         
//                     typename Elements::LocalMatrix Ajac, AjacAnt;
//                     typename Elements::LocalVector Rhs;
//                     typename Elements::Connectivity connec, connecC;
                    
//                     connec = elementsFine_[jel] -> getConnectivity();
                    
//                     //Counts number of coarse mesh intersecting the 
//                     //fine element
//                     int numberIntPoints = elementsFine_[jel] -> 
//                         getNumberOfIntegrationPoints();
//                     int aux;
                    
//                     std::vector<int> ele, diffElem;
//                     ele.reserve(3);
//                     for (int i=0; i<numberIntPoints; i++){
//                         aux = elementsFine_[jel] -> 
//                             getIntegPointCorrespondenceElement(i);
//                         ele.push_back(aux);
//                     };
                    
//                     int numElemIntersect = 1;
//                     int flag = 0;
//                     diffElem.push_back(ele[0]);
                    
//                     for (int i = 1; i<numberIntPoints; i++){
//                         flag = 0;
//                         for (int j = 0; j<numElemIntersect; j++){
//                             if (ele[i] == diffElem[j]) {
//                                 break;
//                             }else{
//                                 flag++;
//                             };
//                             if(flag == numElemIntersect){
//                                 numElemIntersect++;
//                                 diffElem.push_back(ele[i]);
//                             };
//                         };
//                     };
                    
//                     //Compute the Lagrange Multiplier element matrix
//                     for (int ielem = 0; ielem < numElemIntersect; ielem++){
                        
//                         int iElemCoarse = diffElem[ielem];
                        
//                         elementsFine_[jel] -> 
//                             getLagrangeMultipliersDifferentMesh(iElemCoarse);
                        
//                         //Computes element matrix
//                         Ajac = elementsFine_[jel] -> getJacNRMatrix();
                        
//                         typename Elements::LocalVector rhsLagMult,U_;
                        
//                         connecC = elementsCoarse_[iElemCoarse] -> 
//                             getConnectivity();
                        
//                         U_.clear();
//                         for (int i = 0; i < 6; i++){
//                             U_(2*i  ) = nodesCoarse_[connecC(i)] -> 
//                                 getLagrangeMultiplier(0);
//                             U_(2*i+1) = nodesCoarse_[connecC(i)] -> 
//                                 getLagrangeMultiplier(1);
//                         };
                        
//                         noalias(rhsLagMult) = -prod((Ajac),U_);
                        
                        


//                         std::pair<Elements::LocalVector, 
//                                   Elements::LocalMatrix>  lagMult;
                        
//                         lagMult.first.clear();
//                         lagMult.second.clear();

//                         // lagMult = elementsCoarse_[iElemCoarse] -> 
//                         //     getRhsVectorAndBoundaryConditions(Ajac);
//                         lagMult = elementsFine_[jel] -> 
//                             getRhsVectorAndBoundaryConditions(Ajac);
             
//                         //Disperse local contribution into the global matrix
//                         for (int i=0; i<6; i++){
//                             for (int j=0; j<6; j++){
//                                 if (fabs(lagMult.second(2*i  ,2*j  )) >= 1.e-15){
//                                     int dof_i = 2 * numNodesCoarse 
//                                         + numNodesPCoarse + 2 * numNodesFine 
//                                         + numNodesPFine + 2 * connecC(i);
//                                     int dof_j = 2 * numNodesCoarse 
//                                         + numNodesPCoarse + 2 * connec(j);
//                                     ierr = MatSetValues(A,1,&dof_i,1,&dof_j,
//                                                     &lagMult.second(2*i  ,2*j  )
//                                                     ,ADD_VALUES);
                                
//                                     ierr = MatSetValues(A,1,&dof_j,1,&dof_i,
//                                                     &lagMult.second(2*i  ,2*j  )
//                                                     ,ADD_VALUES);
//                                 };
//                                 if (fabs(lagMult.second(2*i+1,2*j  )) >= 1.e-15){
//                                     int dof_i = 2 * numNodesCoarse 
//                                         + numNodesPCoarse + 2 * numNodesFine 
//                                         + numNodesPFine + 2 * connecC(i) + 1;
//                                     int dof_j = 2 * numNodesCoarse 
//                                         + numNodesPCoarse + 2 * connec(j);
//                                     ierr = MatSetValues(A,1,&dof_i,1,&dof_j,
//                                                     &lagMult.second(2*i+1,2*j  )
//                                                     ,ADD_VALUES);
                                
//                                     ierr = MatSetValues(A,1,&dof_j,1,&dof_i,
//                                                     &lagMult.second(2*i+1,2*j  )
//                                                     ,ADD_VALUES);
//                                 };
//                                 if (fabs(lagMult.second(2*i+1,2*j+1)) >= 1.e-15){
//                                     int dof_i = 2 * numNodesCoarse 
//                                         + numNodesPCoarse + 2 * numNodesFine 
//                                         + numNodesPFine + 2 * connecC(i) + 1;
//                                     int dof_j = 2 * numNodesCoarse 
//                                         + numNodesPCoarse + 2 * connec(j) + 1;
//                                     ierr = MatSetValues(A,1,&dof_i,1,&dof_j,
//                                                     &lagMult.second(2*i+1,2*j+1)
//                                                     ,ADD_VALUES);
//                                     ierr = MatSetValues(A,1,&dof_j,1,&dof_i,
//                                                     &lagMult.second(2*i+1,2*j+1)
//                                                     ,ADD_VALUES);
//                                 };
//                                 if (fabs(lagMult.second(2*i  ,2*j+1)) >= 1.e-15){
//                                     int dof_i = 2 * numNodesCoarse 
//                                         + numNodesPCoarse + 2 * numNodesFine 
//                                         + numNodesPFine + 2 * connecC(i);
//                                     int dof_j = 2 * numNodesCoarse 
//                                         + numNodesPCoarse + 2 * connec(j) + 1;
//                                     ierr = MatSetValues(A,1,&dof_i,1,&dof_j,
//                                                     &lagMult.second(2*i  ,2*j+1)
//                                                     ,ADD_VALUES);
                                
//                                     ierr = MatSetValues(A,1,&dof_j,1,&dof_i,
//                                                     &lagMult.second(2*i  ,2*j+1)
//                                                     ,ADD_VALUES);
//                                 };
//                             };
//                             //Rhs vector
//                             if (fabs(lagMult.first(2*i  )) >= 1.e-15){
//                                 int d_i = 2 * numNodesCoarse 
//                                     + numNodesPCoarse + 2 * numNodesFine 
//                                     + numNodesPFine + 2 * connecC(i);
//                                 ierr = VecSetValues(b,1,&d_i,
//                                                     &lagMult.first(2*i  ),
//                                                     ADD_VALUES);
//                             };
//                             if (fabs(lagMult.first(2*i+1)) >= 1.e-15){
//                                 int d_i = 2 * numNodesCoarse 
//                                     + numNodesPCoarse + 2 * numNodesFine 
//                                     + numNodesPFine + 2 * connecC(i) + 1;
//                                 ierr = VecSetValues(b,1,&d_i,
//                                                     &lagMult.first(2*i+1),
//                                                     ADD_VALUES);
//                             };     
                            
                            
//                             if (fabs(rhsLagMult(2*i  )) >= 1.e-15){
//                                 int dof_i = 2 * numNodesCoarse + numNodesPCoarse
//                                     + 2 * connec(i);
//                                 ierr = VecSetValues(b,1,&dof_i,
//                                                     &rhsLagMult(2*i  ),
//                                                     ADD_VALUES);
//                             };
//                             if (fabs(rhsLagMult(2*i+1)) >= 1.e-15){
//                                 int dof_i = 2 * numNodesCoarse + numNodesPCoarse
//                                     + 2 * connec(i) + 1;
//                                 ierr = VecSetValues(b,1,&dof_i,
//                                                     &rhsLagMult(2*i+1),
//                                                     ADD_VALUES);
//                             };
//                         };                                 
//                     }; //Number of intersections
//                 }; // if element belongs to the glue zone
//             }; // Glue zone
//             //------------------------------------------------------------------
//             //-------------------END OF LINEAR SYSTEM ASSEMBLY------------------
//             //------------------------------------------------------------------
            
//             //Assemble matrices and vectors
//             ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
//             ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
            
//             ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
//             ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
            
//             // ierr = MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
//             // ierr = VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
            
//             //Create KSP context to solve the linear system
//             ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
            
//             ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
            



//             // ierr = MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE,0, NULL, &nullsp);
//             // ierr = MatSetNullSpace(A, nullsp);
//             // ierr = MatNullSpaceDestroy(&nullsp);




//             ierr = KSPSetTolerances(ksp,1.e-7,1.e-10,PETSC_DEFAULT,
//                                     100);CHKERRQ(ierr);
            
//             ierr = KSPGMRESSetRestart(ksp, 100); CHKERRQ(ierr);
            
//             ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
            
//             ierr = PCSetType(pc,PCASM);CHKERRQ(ierr);
            
//             //ierr = KSPSetPCSide(ksp, PC_RIGHT);
//             ierr = KSPSetType(ksp,KSPLSQR); CHKERRQ(ierr);
            
//             ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
//             //ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
            



// // #if defined(PETSC_HAVE_MUMPS)
// //             ierr = KSPSetType(ksp,KSPPREONLY);
// //             ierr = KSPGetPC(ksp,&pc);
// //             ierr = PCSetType(pc, PCLU);
// // #endif
            
// //             ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
// //             ierr = KSPSetUp(ksp);
            








//             ierr = KSPSolve(ksp,b,u);CHKERRQ(ierr);
            
//             ierr = KSPGetTotalIterations(ksp, &iterations);
                        
//             //ierr = VecView(u,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
            
//             //Gathers the solution vector to the master process
//             ierr = VecScatterCreateToAll(u, &ctx, &All);CHKERRQ(ierr);
            
//             ierr = VecScatterBegin(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);
//             CHKERRQ(ierr);
            
//             ierr = VecScatterEnd(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);
//             CHKERRQ(ierr);
            
//             ierr = VecScatterDestroy(&ctx);CHKERRQ(ierr);
            
//             //Updates nodal values
//             double u_[2];
//             double norm = 0.;
//             double p_;
//             Ione = 1;
            
//             for (int i = 0; i < numNodesCoarse; ++i){
//                 Ii = 2 * i;
//                 ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
//                 u_[0] = val;
//                 norm += val*val;
//                 nodesCoarse_[i] -> incrementVelocity(0,u_[0]);
                
//                 Ii = 2 * i + 1;
//                 ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
//                 u_[1] = val;
//                 norm += val*val;
//                 nodesCoarse_[i] -> incrementVelocity(1,u_[1]);
//             };
            
//             for (int i = 0; i<numNodesPCoarse; i++){
//                 Ii = 2 * numNodesCoarse + i;
//                 ierr = VecGetValues(All,Ione,&Ii,&val);CHKERRQ(ierr);
//                 p_ = val;
//                 nodesPCoarse_[i] -> incrementPressure(p_);
//             };
            
//             for (int i = 0; i < numNodesFine; ++i){
//                 Ii = 2 * numNodesCoarse + numNodesPCoarse + 2 * i;
//                 ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
//                 u_[0] = val;
//                 norm += val*val;
//                 nodesFine_[i] -> incrementVelocity(0,u_[0]);
                
//                 Ii = 2 * numNodesCoarse + numNodesPCoarse + 2 * i + 1;
//                 ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
//                 u_[1] = val;
//                 norm += val*val;
//                 nodesFine_[i] -> incrementVelocity(1,u_[1]);
//             };
            
//             for (int i = 0; i<numNodesPFine; i++){
//                 Ii = 2* numNodesCoarse + numNodesPCoarse + 2 * numNodesFine + i;
//                 ierr = VecGetValues(All,Ione,&Ii,&val);CHKERRQ(ierr);
//                 p_ = val;
//                 nodesPFine_[i] -> incrementPressure(p_);
//             };
            
//             for (int i = 0; i < numNodesCoarse; ++i){
//                 Ii = 2 * numNodesCoarse + numNodesPCoarse + 2 * numNodesFine 
//                     + numNodesPFine + 2 * i;
//                 ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
//                 u_[0] = val;
//                 nodesCoarse_[i] -> incrementLagrangeMultiplier(0,u_[0]);
                
//                 Ii = 2 * numNodesCoarse + numNodesPCoarse + 2 * numNodesFine 
//                     + numNodesPFine + 2 * i + 1;
//                 ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
//                 u_[1] = val;
//                 nodesCoarse_[i] -> incrementLagrangeMultiplier(1,u_[1]);
                
//                 //std::cout << "LAG M " << u_[0] << " " << u_[1] << std::endl;
//             };
            
//             //Computes the solution vector norm
//             ierr = VecNorm(u,NORM_2,&val);CHKERRQ(ierr);
            
//             boost::posix_time::ptime t2 =
//                 boost::posix_time::microsec_clock::local_time();
            
//             if(rank == 0){
//                 boost::posix_time::time_duration diff = t2 - t1;
                
//                 std::cout<<"Iteration = " << inewton <<  
//                     "      Du Norm = " << std::scientific << sqrt(norm) 
//                          << " " << val << 
//                     "  Time (s) = " << std::fixed << 
//                     diff.total_milliseconds()/1000. << std::endl;
//             };
            
//             ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
//             ierr = VecDestroy(&b); CHKERRQ(ierr);
//             ierr = VecDestroy(&u); CHKERRQ(ierr);
//             ierr = VecDestroy(&All); CHKERRQ(ierr);
//             ierr = MatDestroy(&A); CHKERRQ(ierr);
            
//             if(val <= tolerance){
//                 break;            
//             };          

//             //Updates SUPG Parameter
//             for (int i = 0; i < numElemFine; i++){
//                 elementsFine_[i] -> getParameterSUPG();
//             };   
//             for (int i = 0; i < numElemCoarse; i++){
//                 elementsCoarse_[i] -> getParameterSUPG();
//             };

//         };

//         //Compute real velocity
        
//         QuadShapeFunction<2>                       shapeQuad;
//         typename QuadShapeFunction<2>::Values      phi_;
        
//         for (int i = 0; i<numNodesFine; i++){
//             nodesFine_[i] -> setVelocityArlequin(0,nodesFine_[i] -> 
//                                                  getVelocity(0));
//             nodesFine_[i] -> setVelocityArlequin(1,nodesFine_[i] -> 
//                                                  getVelocity(1));
//         };
        
//         for (int i = 0; i<numNodesGlueZoneFine; i++){
//             typename Nodes::VecLocD xsi;
//             typename Quadrature::NodalValuesQuad u_coarse, v_coarse;
//             typename Elements::Connectivity connecCoarse;
            
//             double u = 0.;
//             double v = 0.;
            
//             int elCoarse = nodesFine_[nodesGlueZoneFine_[i]] -> 
//                 getNodalElemCorrespondence();
//             xsi = nodesFine_[nodesGlueZoneFine_[i]] -> 
//                 getNodalXsiCorrespondence();
            
//             connecCoarse = elementsCoarse_[elCoarse] -> getConnectivity();
            
//             for (int j=0; j<6; j++){
//                 u_coarse(j) = nodesCoarse_[connecCoarse(j)] -> getVelocity(0);
//                 v_coarse(j) = nodesCoarse_[connecCoarse(j)] -> getVelocity(1);
//             };
            
//             shapeQuad.evaluate(xsi,phi_);
            
//             for (int j=0; j<6; j++){
//                 u += u_coarse(j) * phi_(j);
//                 v += v_coarse(j) * phi_(j);
//             };
            
//             double wFunc = nodesFine_[nodesGlueZoneFine_[i]] -> 
//                 getWeightFunction();
            
//             double u_int = nodesFine_[nodesGlueZoneFine_[i]] -> 
//                 getVelocity(0) * wFunc + u * (1. - wFunc);
//             double v_int = nodesFine_[nodesGlueZoneFine_[i]] -> 
//                 getVelocity(1) * wFunc + v * (1. - wFunc);
            
//             nodesFine_[nodesGlueZoneFine_[i]] -> setVelocityArlequin(0,u_int);
//             nodesFine_[nodesGlueZoneFine_[i]] -> setVelocityArlequin(1,v_int);
            
//         };
        
//         for (int i=0; i<numNodesCoarse; i++){
//             nodesCoarse_[i] -> 
//                 setVelocityArlequin(0,nodesCoarse_[i] -> getVelocity(0));
//             nodesCoarse_[i] -> 
//                 setVelocityArlequin(1,nodesCoarse_[i] -> getVelocity(1));
//         };
        
//         //Compute real pressure
        
//         LinShapeFunction<2>                       shapeLin;
//         typename LinShapeFunction<2>::Values      phip_;
        
//         for (int i=0; i<numNodesPFine; i++){
//             nodesPFine_[i] -> 
//                 setPressureArlequin(nodesPFine_[i] -> getPressure());
//         };
        
//         for (int jel = 0; jel<numElemGlueZoneFine; jel++){
            
//             typename Nodes::VecLocD xsi;
//             typename Quadrature::NodalValuesLin p_coarse;
//             typename Elements::ConnectivityP connecPCoarse, connecPFine;
//             typename Elements::Connectivity connecCoarse, connec;
            
//             connec = elementsFine_[elementsGlueZoneFine_[jel]] -> 
//                 getConnectivity();
//             connecPFine = elementsFine_[elementsGlueZoneFine_[jel]] -> 
//                 getConnectivityP();
            
//             for (int i=0; i<3; i++){
                
//                 double p = 0.;
                
//                 int elCoarse = nodesFine_[connec(i)] ->
//                     getNodalElemCorrespondence();
//                 xsi = nodesFine_[connec(i)] -> getNodalXsiCorrespondence();
                
//                 connecCoarse = elementsCoarse_[elCoarse] -> getConnectivity();
//                 connecPCoarse = elementsCoarse_[elCoarse] -> getConnectivityP();
                
//                 for (int j=0; j<3; j++){
//                     p_coarse(j) = nodesPCoarse_[connecPCoarse(j)] -> 
//                         getPressure();
//                 };
                
//                 shapeLin.evaluate(xsi,phip_);
                
//                 for (int j=0; j<3; j++){
//                     p += p_coarse(j) * phip_(j);
//                 };
                
//                 //double wFunc = nodesFine_[connec(i)] -> getWeightFunction();
                
//                 typename Nodes::VecLocD x;

//                 double p_int = 0.;
       
//                 if (fabs(nodesPFine_[connecPFine(i)] -> getPressure()) > 0){
//                     p_int = nodesPFine_[connecPFine(i)] -> getPressure();
//                 }else{
//                     p_int = nodesPFine_[connecPFine(i)] -> getPressure() + p;
//                 };

//                 nodesPFine_[connecPFine(i)] -> setPressureArlequin(p_int);
//             };
//         };
        
//         for (int i=0; i<numNodesPCoarse; i++){
//             nodesPCoarse_[i] -> 
//                 setPressureArlequin(nodesPCoarse_[i] -> getPressure());
//         };
        
        
//         if (rank == 0) {
//             //Computing velocity divergent
//             printVelocity(iTimeStep);
//             printPressure(iTimeStep);
//         };

    // };
        
    return 0;

};


//------------------------------------------------------------------------------
//----------------COMPUTE ARLEQUIN COUPLED NAVIER-STOKES PROBLEM----------------
//------------------------------------------------------------------------------
template<>
int Arlequin<2>::solveArlequinProblemMoving(int iterNumber, double tolerance,
                                            int problem_type, 
                                            int time_dependency){

    Mat               A;
    Vec               b, u, All;
    PetscErrorCode    ierr;
    PetscInt          Istart, Iend, Ii, Ione, iterations;
    KSP               ksp; 
    PC                pc;
    VecScatter        ctx;
    PetscScalar       val;
    //MatNullSpace      nullsp;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);      

    if ((problem_type < 1) || (problem_type > 2)){
        std::cout << "WRONG PROBLEM TYPE." << std::endl;
        return 0;
    };

    if (time_dependency == 0) numTimeSteps = 1;
    
    iTimeStep = 0.;
    
    //Construct the glue zone based on some defined criterion
    setCouplingZone();

    //Computes the Weight function for all the finite elements
    setWeightFunction(4.);
   
    //Computes the Nodal correspondence between fine nodes and coarse elements
    setNodalCorrespondenceFine();

    //Construct the free zone in coarse model
    setFreeZone();

    for (int jel = 0; jel < numElemGlueZoneFine; jel++){
        int iel = elementsGlueZoneFine_[jel];
        //elementsFine_[iel] -> setCompressibility();
    };

    // Computes the system size
    int sysSize = 3 * numNodesCoarse +  
        + 3 * numNodesFine + 2 * numNodesGlueZoneFine;
    
    for (iTimeStep = 0; iTimeStep < numTimeSteps; iTimeStep++){
        
        if (rank == 0) {std::cout << "------------------------- TIME STEP = "
                                  << iTimeStep << " -------------------------"
                                  << std::endl;}
        
        //Updates velocity and acceleration
        for (int i = 0; i < numNodesCoarse; i++){
            double accel[2], u[2], uprev[2];
            
            //Compute acceleration
            u[0] = nodesCoarse_[i] -> getVelocity(0);
            u[1] = nodesCoarse_[i] -> getVelocity(1);
            
            uprev[0] = nodesCoarse_[i] -> getPreviousVelocity(0);
            uprev[1] = nodesCoarse_[i] -> getPreviousVelocity(1);
            
            accel[0] = (u[0] - uprev[0]) / dTime;
            accel[1] = (u[1] - uprev[1]) / dTime;
            
            nodesCoarse_[i] -> setAcceleration(accel);
            
            //Updates velocity
            nodesCoarse_[i] -> setPreviousVelocity(u);
                
        };

        for (int i = 0; i < numNodesFine; i++){
            double accel[2], u[2], uprev[2];
            
            //Compute acceleration
            u[0] = nodesFine_[i] -> getVelocity(0);
            u[1] = nodesFine_[i] -> getVelocity(1);
            
            uprev[0] = nodesFine_[i] -> getPreviousVelocity(0);
            uprev[1] = nodesFine_[i] -> getPreviousVelocity(1);
            
            accel[0] = (u[0] - uprev[0]) / dTime;
            accel[1] = (u[1] - uprev[1]) / dTime;
            
            nodesFine_[i] -> setAcceleration(accel);
            
            //Updates velocity
            nodesFine_[i] -> setPreviousVelocity(u);      

            double T = 50.;
            double h = .5;
            u[0] = 0.;
            u[1] = h * pi / T * cos(pi * iTimeStep * dTime / T);
            nodesFine_[i] -> setMeshVelocity(u);
      
            typename Nodes::VecLocD x, x_ini;
            x = nodesFine_[i] -> getCoordinates();
            x_ini = nodesFine_[i] -> getInitialCoordinates();
            nodesFine_[i] -> setPreviousCoordinates(0,x(0));
            nodesFine_[i] -> setPreviousCoordinates(1,x(1));
            x(1) = x_ini(1) + h * sin(pi * iTimeStep * dTime / T);
            nodesFine_[i] -> setCoordinates(x);
        };

        setFreeZone();
        setSignaledDistance();
        setWeightFunction(4.);
        setNodalCorrespondenceFine();

        //STARTS NEWTON-RAPHSON
        for (int inewton = 0; inewton < iterNumber; inewton++){
            
            boost::posix_time::ptime t1 =
                boost::posix_time::microsec_clock::local_time();
            
            // Preallocates the matrix
            ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                                sysSize, sysSize, 120, NULL, 600, NULL, &A); 
            CHKERRQ(ierr);
            
            // Divides the matrix between the processes
            ierr = MatGetOwnershipRange(A, &Istart, &Iend);CHKERRQ(ierr);
            
            //Create PETSc vectors
            ierr = VecCreate(PETSC_COMM_WORLD, &b); CHKERRQ(ierr);
            ierr = VecSetSizes(b, PETSC_DECIDE, sysSize); CHKERRQ(ierr);
            ierr = VecSetFromOptions(b); CHKERRQ(ierr); 
            ierr = VecDuplicate(b, &u); CHKERRQ(ierr);
            ierr = VecDuplicate(b, &All); CHKERRQ(ierr);
                        
            //------------------------------------------------------------------
            //------------------BEGIN OF LINEAR SYSTEM ASSEMBLY-----------------
            //------------------------------------------------------------------

            //Coarse mesh
            
            for (int jel = 0; jel < numElemCoarse; jel++){   
                
                if (domDecompCoarse.first[jel] == rank) {
                    
                    //Compute Element matrix
                    if (time_dependency == 0){ 
                        if (problem_type == 1){
                            elementsCoarse_[jel] -> getSteadyStokes();
                        }else{
                            if (inewton == 0){
                                elementsCoarse_[jel] -> getSteadyStokes();
                            }else{
                                elementsCoarse_[jel] -> getSteadyNavierStokes();
                            };
                        };
                    }else{
                        if (problem_type == 1){
                            elementsCoarse_[jel] -> getTransientNavierStokes();
                        }else{
                            elementsCoarse_[jel] -> getTransientNavierStokes();
                        };
                    };

                    typename Elements::LocalMatrix Ajac;
                    typename Elements::LocalVector Rhs;
                    typename Elements::Connectivity connec;
                    
                    //Gets element connectivity, jacobian and rhs 
                    connec = elementsCoarse_[jel] -> getConnectivity();
                    Ajac = elementsCoarse_[jel] -> getJacNRMatrix();
                    Rhs = elementsCoarse_[jel] -> getRhsVector();
                    
                    //Disperse local contributions into the global matrix
                    for (int i=0; i<6; i++){
                        for (int j=0; j<6; j++){
                            
                            if (fabs(Ajac(2*i  ,2*j  )) >= 1.e-15){
                                int dof_i = 2 * connec(i);
                                int dof_j = 2 * connec(j);
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    
                                                    &Ajac(2*i  ,2*j  ),
                                                    ADD_VALUES);
                            };
                            if (fabs(Ajac(2*i+1,2*j  )) >= 1.e-15){
                                int dof_i = 2 * connec(i) + 1;
                                int dof_j = 2 * connec(j);
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    
                                                    &Ajac(2*i+1,2*j  ),
                                                    ADD_VALUES);
                            };

                            if (fabs(Ajac(2*i  ,2*j+1)) >= 1.e-15){
                                int dof_i = 2 * connec(i);
                                int dof_j = 2 * connec(j) + 1;
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    
                                                    &Ajac(2*i  ,2*j+1),
                                                    ADD_VALUES);
                            };

                            if (fabs(Ajac(2*i+1,2*j+1)) >= 1.e-15){
                                int dof_i = 2 * connec(i) + 1;
                                int dof_j = 2 * connec(j) + 1;
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    
                                                    &Ajac(2*i+1,2*j+1),
                                                    ADD_VALUES);
                            };

                            //Matrix Q and Qt
                            if (fabs(Ajac(2*i  ,12+j)) >= 1.e-15){
                                int dof_i = 2 * connec(i);
                                int dof_j = 2 * numNodesCoarse + connec(j);
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,
                                                    &Ajac(2*i  ,12+j),
                                                    ADD_VALUES);
                            };
                            
                            if (fabs(Ajac(12+j,2*i  )) >= 1.e-15){
                                int dof_i = 2 * connec(i);
                                int dof_j = 2 * numNodesCoarse + connec(j);
                                ierr = MatSetValues(A,1,&dof_j,1,&dof_i,    
                                                    &Ajac(12+j,2*i  ),
                                                    ADD_VALUES);
                            };

                            if (fabs(Ajac(2*i+1,12+j)) >= 1.e-15){
                                int dof_i = 2 * connec(i) + 1;
                                int dof_j = 2 * numNodesCoarse + connec(j);
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,
                                                    &Ajac(2*i+1,12+j),
                                                    ADD_VALUES);
                            };

                            if (fabs(Ajac(12+j,2*i+1)) >= 1.e-15){
                                int dof_i = 2 * connec(i) + 1;
                                int dof_j = 2 * numNodesCoarse + connec(j);
                                ierr = MatSetValues(A,1,&dof_j,1,&dof_i,    
                                                    &Ajac(12+j,2*i+1),
                                                    ADD_VALUES);
                            };

                            if (fabs(Ajac(12+i,12+j)) >= 1.e-15){
                                int dof_i = 2 * numNodesCoarse + connec(i);
                                int dof_j = 2 * numNodesCoarse + connec(j);
                                ierr = MatSetValues(A,1,&dof_j,1,&dof_i,    
                                                    &Ajac(12+i,12+j),
                                                    ADD_VALUES);
                            };
                        };                  
                        //Rhs vector
                        if (fabs(Rhs(2*i  )) >= 1.e-15){
                            int dof_i = 2 * connec(i);
                            ierr = VecSetValues(b,1,&dof_i,&Rhs(2*i  ),
                                                ADD_VALUES);
                        };

                        if (fabs(Rhs(2*i+1)) >= 1.e-15){
                            int dof_i = 2 * connec(i) + 1;
                            ierr = VecSetValues(b,1,&dof_i,&Rhs(2*i+1),
                                                ADD_VALUES);
                        };

                        if (fabs(Rhs(12+i)) >= 1.e-15){
                            int dof_i = 2 * numNodesCoarse + connec(i);
                            ierr = VecSetValues(b,1,&dof_i,&Rhs(12+i),
                                                ADD_VALUES);
                        };
                    };
                };
            };

            //Fine mesh
            for (int jel = 0; jel < numElemFine; jel++){   
                
                if (domDecompFine.first[jel] == rank) {
                    
                    //Compute Element matrix
                    if (time_dependency == 0) {
                        if (problem_type == 1){
                            elementsFine_[jel] -> getSteadyStokes();
                        }else{
                            if (inewton == 0){
                                elementsFine_[jel] -> getSteadyStokes();
                            }else{
                                elementsFine_[jel] -> getSteadyNavierStokes();
                            };
                        };
                    }else{
                        if (problem_type == 1){
                            elementsFine_[jel] -> getTransientNavierStokes();
                        }else{
                            elementsFine_[jel] -> getTransientNavierStokes();
                        };
                    };

                    typename Elements::LocalMatrix Ajac;
                    typename Elements::LocalVector Rhs;
                    typename Elements::Connectivity connec;
                    
                    //Gets element connectivity, jacobian and rhs 
                    connec = elementsFine_[jel] -> getConnectivity();
                    Ajac = elementsFine_[jel] -> getJacNRMatrix();
                    Rhs = elementsFine_[jel] -> getRhsVector();              

                    //Disperse local contributions into the global matrix
                    for (int i=0; i<6; i++){
                        for (int j=0; j<6; j++){
                            if (fabs(Ajac(2*i  ,2*j  )) >= 1.e-15){
                                int dof_i = 3 * numNodesCoarse + 2 * connec(i);
                                int dof_j = 3 * numNodesCoarse + 2 * connec(j);
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    
                                                    &Ajac(2*i  ,2*j  ),
                                                    ADD_VALUES);
                            };
                            
                            if (fabs(Ajac(2*i+1,2*j  )) >= 1.e-15){
                                int dof_i = 3 * numNodesCoarse +
                                    2 * connec(i) + 1;
                                int dof_j = 3 * numNodesCoarse + 2 * connec(j);
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    
                                                    &Ajac(2*i+1,2*j  ),
                                                    ADD_VALUES);
                            };

                            if (fabs(Ajac(2*i  ,2*j+1)) >= 1.e-15){
                                int dof_i = 3 * numNodesCoarse + 2 * connec(i);
                                int dof_j = 3 * numNodesCoarse + 
                                    2 * connec(j) + 1;
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    
                                                    &Ajac(2*i  ,2*j+1),
                                                    ADD_VALUES);
                            };

                            if (fabs(Ajac(2*i+1,2*j+1)) >= 1.e-15){
                                int dof_i = 3 * numNodesCoarse + 
                                    2 * connec(i) + 1;
                                int dof_j = 3 * numNodesCoarse + 
                                    2 * connec(j) + 1;
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    
                                                    &Ajac(2*i+1,2*j+1),
                                                    ADD_VALUES);
                            };

                            //Matrix Q and Qt
                            if (fabs(Ajac(2*i  ,12+j)) >= 1.e-15){
                                int dof_i = 3 * numNodesCoarse + 2 * connec(i);
                                int dof_j = 3 * numNodesCoarse + 
                                    2 * numNodesFine + connec(j);
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,
                                                    &Ajac(2*i  ,12+j),
                                                    ADD_VALUES);
                            };

                            if (fabs(Ajac(12+j,2*i  )) >= 1.e-15){
                                int dof_i = 3 * numNodesCoarse + 2 * connec(i);
                                int dof_j = 3 * numNodesCoarse + 
                                    2 * numNodesFine + connec(j);
                                ierr = MatSetValues(A,1,&dof_j,1,&dof_i,  
                                                    &Ajac(12+j,2*i  ),
                                                    ADD_VALUES);
                            };

                            if (fabs(Ajac(2*i+1,12+j)) >= 1.e-15){
                                int dof_i = 3 * numNodesCoarse + 
                                    2 * connec(i) + 1;
                                int dof_j = 3 * numNodesCoarse + 
                                    2 * numNodesFine + connec(j);
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,
                                                    &Ajac(2*i+1,12+j),
                                                    ADD_VALUES);
                            };
                            
                            if (fabs(Ajac(12+j,2*i+1)) >= 1.e-15){
                                int dof_i = 3 * numNodesCoarse + 
                                    2 * connec(i) + 1;
                                int dof_j = 3 * numNodesCoarse + 
                                    2 * numNodesFine + connec(j);
                                ierr = MatSetValues(A,1,&dof_j,1,&dof_i,  
                                                    &Ajac(12+j,2*i+1),
                                                    ADD_VALUES);
                            };
                            
                            if (fabs(Ajac(12+i,12+j)) >= 1.e-15){
                                int dof_i = 3 * numNodesCoarse + 
                                    2 * numNodesFine + connec(i);
                                int dof_j = 3 * numNodesCoarse + 
                                    2 * numNodesFine + connec(j);
                                ierr = MatSetValues(A,1,&dof_j,1,&dof_i,  
                                                    &Ajac(12+i,12+j),
                                                    ADD_VALUES);
                            };
                        };     
                        ///Rhs vector
                        if (fabs(Rhs(2*i  )) >= 1.e-15){
                            int dof_i = 3 * numNodesCoarse + 2 * connec(i);
                            ierr = VecSetValues(b,1,&dof_i,&Rhs(2*i  ),
                                                ADD_VALUES);
                        };

                        if (fabs(Rhs(2*i+1)) >= 1.e-15){
                            int dof_i = 3 * numNodesCoarse + 2 * connec(i) + 1;
                            ierr = VecSetValues(b,1,&dof_i,&Rhs(2*i+1),
                                                ADD_VALUES);
                        };

                        if (fabs(Rhs(12+i)) >= 1.e-15){
                            int dof_i = 3 * numNodesCoarse + 2 * numNodesFine +
                                connec(i);
                            ierr = VecSetValues(b,1,&dof_i,&Rhs(12+i),
                                                ADD_VALUES);
                        };
                    }; 
                };                
            };
              

            //Lagrange Multipliers
            for (int l=0; l< numElemGlueZoneFine; l++){
                
                int jel = elementsGlueZoneFine_[l];
                
                if (domDecompFine.first[jel] == rank) {
                    
                    typename Elements::LocalMatrix Ajac, AjacAnt, AStab;
                    typename Elements::LocalVector Rhs, RhsStab;
                    typename Elements::Connectivity connec, connecC, connecL;
                    
                    connec = elementsFine_[jel] -> getConnectivity();
                    connecL = glueZoneFine_[l] -> getConnectivity();

                    // FINE MESH
                    
                    //Computes element matrix
                    elementsFine_[jel] -> getLagrangeMultipliersSameMesh();
                    
                    //Gets element matrice and rhs vectors
                    Ajac = - elementsFine_[jel] -> getLagrMultMatrix();
                    Rhs = - elementsFine_[jel] -> getRhsVector();
                    // AjacAnt = - elementsFine_[jel] -> getJacNRMatrix();
                    Ajac = trans(Ajac);

                    // elementsFine_[jel] -> getLMStabilizationSameMesh();

                    // AStab = - elementsFine_[jel] -> getJacNRMatrix();
                    // RhsStab = - elementsFine_[jel] -> getRhsVector();

                    //AjacAnt = trans(AjacAnt);
                    
                    typename Elements::LocalVector rhsLagMult,U_,lagStab;
                    
                    U_.clear();
                    rhsLagMult.clear();
                    lagStab.clear();
                    for (int i = 0; i < 6; i++){
                        U_(2*i  ) = nodesFine_[connec(i)] -> 
                            getLagrangeMultiplier(0);
                        U_(2*i+1) = nodesFine_[connec(i)] -> 
                            getLagrangeMultiplier(1);
                    };
                    
                    noalias(rhsLagMult) = -prod(trans(Ajac),U_);
                    // noalias(lagStab) = -prod(trans(AjacAnt),U_);
                   
                    
                    //Disperse local contributions into the global matrix
                    for (int i=0; i<6; i++){
                        for (int j=0; j<6; j++){
                            
                            if (fabs(Ajac(2*i  ,2*j  )) >= 1.e-15){
                                int d_i = 3 * numNodesCoarse + 3 * numNodesFine
                                    + 2 * connecL(i);
                                int d_j = 3 * numNodesCoarse + 2 * connec(j);
                                ierr = MatSetValues(A,1,&d_i,1,&d_j,
                                                    &Ajac(2*i  ,2*j  ),
                                                    ADD_VALUES);
                                ierr = MatSetValues(A,1,&d_j,1,&d_i,
                                                    &Ajac(2*i  ,2*j  ),
                                                    ADD_VALUES);
                            };
                            // ierr = MatSetValues(A,1,&d_i,1,&d_j,
                            //           &AStab(2*i  ,2*j  ),ADD_VALUES);

                            if (fabs(Ajac(2*i+1,2*j  )) >= 1.e-15){
                                int d_i = 3 * numNodesCoarse + 3 * numNodesFine
                                    + 2 * connecL(i) + 1;
                                int d_j = 3 * numNodesCoarse + 2 * connec(j);
                                ierr = MatSetValues(A,1,&d_i,1,&d_j,
                                                    &Ajac(2*i+1,2*j  ),
                                                    ADD_VALUES);
                                ierr = MatSetValues(A,1,&d_j,1,&d_i,
                                                    &Ajac(2*i+1,2*j  ),
                                                    ADD_VALUES);
                            };
                            // ierr = MatSetValues(A,1,&d_i,1,&d_j,
                            //                     &AStab(2*i+1,2*j  ),ADD_VALUES);

                            if (fabs(Ajac(2*i+1,2*j+1)) >= 1.e-15){
                                int d_i = 3 * numNodesCoarse + 3 * numNodesFine
                                    + 2 * connecL(i) + 1;
                                int d_j = 3 * numNodesCoarse + 2*connec(j) + 1;
                                ierr = MatSetValues(A,1,&d_i,1,&d_j,
                                                    &Ajac(2*i+1,2*j+1),
                                                    ADD_VALUES);
                                ierr = MatSetValues(A,1,&d_j,1,&d_i,
                                                    &Ajac(2*i+1,2*j+1),
                                                    ADD_VALUES);
                            };
                            // ierr = MatSetValues(A,1,&d_i,1,&d_j,
                            //                     &AStab(2*i+1,2*j+1),ADD_VALUES);

                            if (fabs(Ajac(2*i  ,2*j+1)) >= 1.e-15){
                                int d_i = 3 * numNodesCoarse + 3 * numNodesFine
                                    + 2 * connecL(i);
                                int d_j = 3 * numNodesCoarse + 2*connec(j) + 1;
                                ierr = MatSetValues(A,1,&d_i,1,&d_j,
                                                    &Ajac(2*i  ,2*j+1),
                                                    ADD_VALUES);
                                ierr = MatSetValues(A,1,&d_j,1,&d_i,
                                                    &Ajac(2*i  ,2*j+1),
                                                    ADD_VALUES);
                            };
                            // ierr = MatSetValues(A,1,&d_i,1,&d_j,
                            //                     &AStab(2*i  ,2*j+1),ADD_VALUES);
                            // if (fabs(AjacAnt(2*i  ,2*j  )) >= 1.e-15){
                            //     int d_i = 3 * numNodesCoarse + 3 * numNodesFine
                            //         + 2 * connecL(i);
                            //     int d_j = 3 * numNodesCoarse + 3 * numNodesFine
                            //         + 2 * connecL(j);
                            //     ierr = MatSetValues(A,1,&d_i,1,&d_j,
                            //                         &AjacAnt(2*i  ,2*j  )
                            //                         ,ADD_VALUES);
                            // };

                            // if (fabs(AjacAnt(2*i+1,2*j+1)) >= 1.e-15){
                            //     int d_i = 3 * numNodesCoarse + 3 * numNodesFine
                            //         + 2 * connecL(i) + 1;
                            //     int d_j = 3 * numNodesCoarse + 3 * numNodesFine
                            //         + 2 * connecL(j) + 1;
                            //     ierr = MatSetValues(A,1,&d_i,1,&d_j,
                            //                         &AjacAnt(2*i+1,2*j+1)
                            //                         ,ADD_VALUES);
                            // };
                        };
                        //Rhs vector

                        int dof_i = 3 * numNodesCoarse + 3 * numNodesFine +
                            2 * connecL(i);
                        ierr = VecSetValues(b,1,&dof_i,&Rhs(2*i  ),
                                            ADD_VALUES);
                        // ierr = VecSetValues(b,1,&dof_i,&RhsStab(2*i  ),
                        //                     ADD_VALUES);
                        // ierr = VecSetValues(b,1,&dof_i,&lagStab(2*i  ),
                        //                     ADD_VALUES);

                        dof_i = 3 * numNodesCoarse + 3 * numNodesFine + 
                            2 * connecL(i) + 1;
                        ierr = VecSetValues(b,1,&dof_i,&Rhs(2*i+1),
                                            ADD_VALUES);
                        // ierr = VecSetValues(b,1,&dof_i,&RhsStab(2*i+1),
                        //                     ADD_VALUES);
                        // ierr = VecSetValues(b,1,&dof_i,&lagStab(2*i+1),
                        //                     ADD_VALUES);

                        dof_i = 3 * numNodesCoarse + 2 * connec(i);
                        ierr = VecSetValues(b,1,&dof_i,&rhsLagMult(2*i  )
                                            ,ADD_VALUES);

                        dof_i = 3 * numNodesCoarse + 2 * connec(i) + 1;
                        ierr = VecSetValues(b,1,&dof_i,&rhsLagMult(2*i+1)
                                            ,ADD_VALUES);
                    };      

                    //COAESE MESH
                    
                    //Counts number of coarse mesh intersecting the 
                    //fine element
                    int numberIntPoints = elementsFine_[jel] -> 
                        getNumberOfIntegrationPoints();
                    int aux;
                    
                    std::vector<int> ele, diffElem;
                    ele.clear();
                    diffElem.clear();
                    ele.reserve(3);
                    for (int i=0; i<numberIntPoints; i++){
                        aux = elementsFine_[jel] -> 
                            getIntegPointCorrespondenceElement(i);
                        ele.push_back(aux);

                    };
                    
                    int numElemIntersect = 1;
                    int flag = 0;
                    diffElem.push_back(ele[0]);
                    
                    for (int i = 1; i<numberIntPoints; i++){
                        flag = 0;
                        for (int j = 0; j<numElemIntersect; j++){
                            if (ele[i] == diffElem[j]) {
                                break;
                            }else{
                                flag++;
                            };
                            if(flag == numElemIntersect){
                                numElemIntersect++;
                                diffElem.push_back(ele[i]);
                            };
                        };
                    };
                    //std::cout << "JEL " << jel << " " << numElemIntersect << std::endl;
                    //Compute the Lagrange Multiplier element matrix
                    for (int ielem = 0; ielem < numElemIntersect; ielem++){
                        
                        int iElemCoarse = diffElem[ielem];
                       
                        elementsFine_[jel] -> 
                            getLagrangeMultipliersDifferentMesh(iElemCoarse);
                        
                        //Computes element matrix
                        Ajac = elementsFine_[jel] -> getLagrMultMatrix();
                        Ajac = trans(Ajac);
                        // AjacAnt = elementsFine_[jel] -> getJacNRMatrix();
                        // AjacAnt = trans(AjacAnt);

                        connecC = elementsCoarse_[iElemCoarse] -> 
                            getConnectivity();

                        
                        typename Elements::LocalVector rhsLagMult,U_,lagStab;
                        
                        U_.clear();
                        for (int i = 0; i < 6; i++){
                            U_(2*i  ) = nodesFine_[connec(i)] -> 
                                getLagrangeMultiplier(0);
                            U_(2*i+1) = nodesFine_[connec(i)] -> 
                                getLagrangeMultiplier(1);
                        };
                        
                        noalias(rhsLagMult) = -prod(trans(Ajac),U_);
                        // noalias(lagStab) = -prod(trans(AjacAnt),U_);
                        
                        std::pair<Elements::LocalVector, 
                                  Elements::LocalMatrix>  lagMult;
                        
                        lagMult.first.clear();
                        lagMult.second.clear();
                         
                        lagMult = elementsCoarse_[iElemCoarse] -> 
                            getRhsVectorAndBoundaryConditions(Ajac);
                        
                        //Disperse local contribution into the global matrix
                        for (int i=0; i<6; i++){
                            for (int j=0; j<6; j++){

                                if (fabs(Ajac(2*i  ,2*j  )) >= 1.e-15){
                                    int dof_i = 3 * numNodesCoarse 
                                        + 3 * numNodesFine + 2 * connecL(i);
                                    int dof_j = 2 * connecC(j);
                                    ierr = MatSetValues(A,1,&dof_i,1,&dof_j,
                                                    &lagMult.second(2*i  ,2*j  )
                                                        ,ADD_VALUES);
                                    ierr = MatSetValues(A,1,&dof_j,1,&dof_i,
                                                    &lagMult.second(2*i  ,2*j  )
                                                        ,ADD_VALUES);
                                };

                                if (fabs(Ajac(2*i+1,2*j  )) >= 1.e-15){
                                    int dof_i = 3 * numNodesCoarse 
                                        + 3 * numNodesFine + 2 * connecL(i) + 1;
                                    int dof_j = 2 * connecC(j);
                                    ierr = MatSetValues(A,1,&dof_i,1,&dof_j,
                                                    &lagMult.second(2*i+1,2*j  )
                                                        ,ADD_VALUES);
                                    ierr = MatSetValues(A,1,&dof_j,1,&dof_i,
                                                    &lagMult.second(2*i+1,2*j  )
                                                        ,ADD_VALUES);
                                };

                                if (fabs(Ajac(2*i+1,2*j+1)) >= 1.e-15){
                                    int dof_i = 3 * numNodesCoarse 
                                        + 3 * numNodesFine + 2 * connecL(i) + 1;
                                    int dof_j = 2 * connecC(j) + 1;
                                    ierr = MatSetValues(A,1,&dof_i,1,&dof_j,
                                                    &lagMult.second(2*i+1,2*j+1)
                                                        ,ADD_VALUES);
                                    ierr = MatSetValues(A,1,&dof_j,1,&dof_i,
                                                    &lagMult.second(2*i+1,2*j+1)
                                                        ,ADD_VALUES);
                                };

                                if (fabs(Ajac(2*i  ,2*j+1)) >= 1.e-15){
                                    int dof_i = 3 * numNodesCoarse 
                                        + 3 * numNodesFine + 2 * connecL(i);
                                    int dof_j = 2 * connecC(j) + 1;
                                    ierr = MatSetValues(A,1,&dof_i,1,&dof_j,
                                                    &lagMult.second(2*i  ,2*j+1)
                                                        ,ADD_VALUES);
                                    ierr = MatSetValues(A,1,&dof_j,1,&dof_i,
                                                    &lagMult.second(2*i  ,2*j+1)
                                                        ,ADD_VALUES);
                                };

                                // if (fabs(AjacAnt(2*i  ,2*j  )) >= 1.e-15){
                                //     int d_i = 3 * numNodesCoarse + 
                                //         3 * numNodesFine + 2 * connecL(i);
                                //     int d_j = 3 * numNodesCoarse + 
                                //         3 * numNodesFine + 2 * connecL(j);
                                //     ierr = MatSetValues(A,1,&d_i,1,&d_j,
                                //                         &AjacAnt(2*i  ,2*j  ),
                                //                         ADD_VALUES);
                                // };
                                // if (fabs(AjacAnt(2*i+1,2*j+1)) >= 1.e-15){
                                //     int d_i = 3 * numNodesCoarse + 
                                //         3 * numNodesFine + 2 * connecL(i) + 1;
                                //     int d_j = 3 * numNodesCoarse + 
                                //         3 * numNodesFine + 2 * connecL(j) + 1;
                                //     ierr = MatSetValues(A,1,&d_i,1,&d_j,
                                //                         &AjacAnt(2*i+1,2*j+1),
                                //                         ADD_VALUES);
                                // };
                            };
                            //Rhs vector

                            int d_i = 3 * numNodesCoarse 
                                + 3 * numNodesFine + 2 * connecL(i);
                            ierr = VecSetValues(b,1,&d_i,
                                                &lagMult.first(2*i  ),
                                                ADD_VALUES);
                            // ierr = VecSetValues(b,1,&d_i,
                            //                     &lagStab(2*i  ),
                            //                     ADD_VALUES);

                            d_i = 3 * numNodesCoarse 
                                + 3 * numNodesFine + 2 * connecL(i) + 1;
                            ierr = VecSetValues(b,1,&d_i,
                                                &lagMult.first(2*i+1),
                                                ADD_VALUES);
                            // ierr = VecSetValues(b,1,&d_i,
                            //                     &lagStab(2*i+1),
                            //                     ADD_VALUES);

                            int dof_i = 2 * connecC(i);
                            ierr = VecSetValues(b,1,&dof_i,
                                                &rhsLagMult(2*i  ),
                                                ADD_VALUES);

                            dof_i = 2 * connecC(i) + 1;
                            ierr = VecSetValues(b,1,&dof_i,
                                                &rhsLagMult(2*i+1),
                                                ADD_VALUES);
                        };                                 
                    }; //Number of intersections
                }; // if element belongs to the glue zone
            }; // Glue zone


            //------------------------------------------------------------------
            //-------------------END OF LINEAR SYSTEM ASSEMBLY------------------
            //------------------------------------------------------------------
            //std::cout << "Enter PETSc " << rank << std::endl;
            
            //Assemble matrices and vectors
            ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
            ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
            
            ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
            ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
            
            //ierr = MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
            // ierr = VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
            
            //Create KSP context to solve the linear system
            ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
            
            ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
            



            // ierr = MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE,0, NULL, &nullsp);
            // ierr = MatSetNullSpace(A, nullsp);
            // ierr = MatNullSpaceDestroy(&nullsp);




            ierr = KSPSetTolerances(ksp,1.e-7,1.e-10,PETSC_DEFAULT,
                                    1000);CHKERRQ(ierr);
            
            //ierr = KSPGMRESSetRestart(ksp, 10); CHKERRQ(ierr);
            
            ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
            
            ierr = PCSetType(pc,PCNONE);CHKERRQ(ierr);
            
            //ierr = KSPSetPCSide(ksp, PC_RIGHT);
            ierr = KSPSetType(ksp,KSPLSQR); CHKERRQ(ierr);
            
            ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
            // ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);



            
// #if defined(PETSC_HAVE_MUMPS)
//             ierr = KSPSetType(ksp,KSPPREONLY);
//             ierr = KSPGetPC(ksp,&pc);
//             ierr = PCSetType(pc, PCLU);
// #endif
            
//             ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
//             ierr = KSPSetUp(ksp);
            




           // ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);

            ierr = KSPSolve(ksp,b,u);CHKERRQ(ierr);
            
            ierr = KSPGetTotalIterations(ksp, &iterations);
            
            //if (rank == 0)std::cout << "GMRES Iterations = " << iterations << std::endl;
            
            //ierr = VecView(u,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
            
            //Gathers the solution vector to the master process
            ierr = VecScatterCreateToAll(u, &ctx, &All);CHKERRQ(ierr);
            
            ierr = VecScatterBegin(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);
            CHKERRQ(ierr);
            
            ierr = VecScatterEnd(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);
            CHKERRQ(ierr);
            
            ierr = VecScatterDestroy(&ctx);CHKERRQ(ierr);
            
            //Updates nodal values
            double u_[2];
            double norm = 0.;
            double p_;
            Ione = 1;
            
            for (int i = 0; i < numNodesCoarse; ++i){
                Ii = 2 * i;
                ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                u_[0] = val;
                norm += val*val;
                nodesCoarse_[i] -> incrementVelocity(0,u_[0]);
                
                Ii = 2 * i + 1;
                ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                u_[1] = val;
                norm += val*val;
                nodesCoarse_[i] -> incrementVelocity(1,u_[1]);
            };
            
            for (int i = 0; i<numNodesCoarse; i++){
                Ii = 2 * numNodesCoarse + i;
                ierr = VecGetValues(All,Ione,&Ii,&val);CHKERRQ(ierr);
                p_ = val;
                nodesCoarse_[i] -> incrementPressure(p_);
            };
            
            for (int i = 0; i < numNodesFine; ++i){
                Ii = 3 * numNodesCoarse + 2 * i;
                ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                u_[0] = val;
                norm += val*val;
                nodesFine_[i] -> incrementVelocity(0,u_[0]);
                
                Ii = 3 * numNodesCoarse + 2 * i + 1;
                ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                u_[1] = val;
                norm += val*val;
                nodesFine_[i] -> incrementVelocity(1,u_[1]);
            };
            
            for (int i = 0; i<numNodesFine; i++){
                Ii = 3 * numNodesCoarse + 2 * numNodesFine + i;
                ierr = VecGetValues(All,Ione,&Ii,&val);CHKERRQ(ierr);
                p_ = val;
                nodesFine_[i] -> incrementPressure(p_);
            };
            
            for (int i = 0; i < numNodesGlueZoneFine; ++i){
                Ii = 3 * numNodesCoarse + 3 * numNodesFine + 2 * i;
                ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                u_[0] = val;
                nodesFine_[nodesGlueZoneFine_[i]] -> incrementLagrangeMultiplier(0,u_[0]);
                
                Ii = 3 * numNodesCoarse + 3 * numNodesFine + 2 * i + 1;
                ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                u_[1] = val;
                nodesFine_[nodesGlueZoneFine_[i]] -> incrementLagrangeMultiplier(1,u_[1]);
                
                // std::cout << "LAG M " << u_[0] << " " << u_[1] << std::endl;
            };
            
            //Computes the solution vector norm
            ierr = VecNorm(u,NORM_2,&val);CHKERRQ(ierr);
            
            boost::posix_time::ptime t2 =
                boost::posix_time::microsec_clock::local_time();
            
            if(rank == 0){
                boost::posix_time::time_duration diff = t2 - t1;
                
                std::cout<<"Iteration = " << inewton << " (" << iterations <<  
                    ")  Du Norm = " << std::scientific << sqrt(norm) 
                         << " " << val << 
                    "  Time (s) = " << std::fixed << 
                    diff.total_milliseconds()/1000. << std::endl;
            };
            
            ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
            ierr = VecDestroy(&b); CHKERRQ(ierr);
            ierr = VecDestroy(&u); CHKERRQ(ierr);
            ierr = VecDestroy(&All); CHKERRQ(ierr);
            ierr = MatDestroy(&A); CHKERRQ(ierr);
            
            if(val <= tolerance){
                break;            
            };          

            //Updates SUPG Parameter
            // for (int i = 0; i < numElemFine; i++){
            //     elementsFine_[i] -> getParameterSUPG();
            // };   
            // for (int i = 0; i < numElemCoarse; i++){
            //     elementsCoarse_[i] -> getParameterSUPG();
            // };

        };

        //Compute real velocity
        
        QuadShapeFunction<2>                       shapeQuad;
        typename QuadShapeFunction<2>::Values      phi_;
        
        for (int i = 0; i<numNodesFine; i++){
            nodesFine_[i] -> setVelocityArlequin(0,nodesFine_[i] -> 
                                                 getVelocity(0));
            nodesFine_[i] -> setVelocityArlequin(1,nodesFine_[i] -> 
                                                 getVelocity(1));
        };
        
        for (int i = 0; i<numNodesGlueZoneFine; i++){
            typename Nodes::VecLocD xsi;
            typename Quadrature::NodalValuesQuad u_coarse, v_coarse;
            typename Elements::Connectivity connecCoarse;
            
            double u = 0.;
            double v = 0.;
            
            int elCoarse = nodesFine_[nodesGlueZoneFine_[i]] -> 
                getNodalElemCorrespondence();
            xsi = nodesFine_[nodesGlueZoneFine_[i]] -> 
                getNodalXsiCorrespondence();
            
            connecCoarse = elementsCoarse_[elCoarse] -> getConnectivity();
            
            for (int j=0; j<6; j++){
                u_coarse(j) = nodesCoarse_[connecCoarse(j)] -> getVelocity(0);
                v_coarse(j) = nodesCoarse_[connecCoarse(j)] -> getVelocity(1);
            };
            
            shapeQuad.evaluate(xsi,phi_);
            
            for (int j=0; j<6; j++){
                u += u_coarse(j) * phi_(j);
                v += v_coarse(j) * phi_(j);
            };
            
            double wFunc = nodesFine_[nodesGlueZoneFine_[i]] -> 
                getWeightFunction();
            
            double u_int = nodesFine_[nodesGlueZoneFine_[i]] -> 
                getVelocity(0) * wFunc + u * (1. - wFunc);
            double v_int = nodesFine_[nodesGlueZoneFine_[i]] -> 
                getVelocity(1) * wFunc + v * (1. - wFunc);
            
            nodesFine_[nodesGlueZoneFine_[i]] -> setVelocityArlequin(0,u_int);
            nodesFine_[nodesGlueZoneFine_[i]] -> setVelocityArlequin(1,v_int);
            
        };
        
        for (int i=0; i<numNodesCoarse; i++){
            nodesCoarse_[i] -> 
                setVelocityArlequin(0,nodesCoarse_[i] -> getVelocity(0));
            nodesCoarse_[i] -> 
                setVelocityArlequin(1,nodesCoarse_[i] -> getVelocity(1));
        };
        
        //Compute real pressure
        
        // LinShapeFunction<2>                       shapeLin;
        // typename LinShapeFunction<2>::Values      phip_;
        
        // for (int i=0; i<numNodesPFine; i++){
        //     nodesPFine_[i] -> 
        //         setPressureArlequin(nodesPFine_[i] -> getPressure());
        // };
        
        // for (int jel = 0; jel<numElemGlueZoneFine; jel++){
            
        //     typename Nodes::VecLocD xsi;
        //     typename Quadrature::NodalValuesLin p_coarse;
        //     typename Elements::ConnectivityP connecPCoarse, connecPFine;
        //     typename Elements::Connectivity connecCoarse, connec;
            
        //     connec = elementsFine_[elementsGlueZoneFine_[jel]] -> 
        //         getConnectivity();
        //     connecPFine = elementsFine_[elementsGlueZoneFine_[jel]] -> 
        //         getConnectivityP();
            
        //     for (int i=0; i<3; i++){
                
        //         double p = 0.;
                
        //         int elCoarse = nodesFine_[connec(i)] ->
        //             getNodalElemCorrespondence();
        //         xsi = nodesFine_[connec(i)] -> getNodalXsiCorrespondence();
                
        //         connecCoarse = elementsCoarse_[elCoarse] -> getConnectivity();
        //         connecPCoarse = elementsCoarse_[elCoarse] -> getConnectivityP();
                
        //         for (int j=0; j<3; j++){
        //             p_coarse(j) = nodesPCoarse_[connecPCoarse(j)] -> 
        //                 getPressure();
        //         };
                
        //         shapeLin.evaluate(xsi,phip_);
                
        //         for (int j=0; j<3; j++){
        //             p += p_coarse(j) * phip_(j);
        //         };
                
        //         //double wFunc = nodesFine_[connec(i)] -> getWeightFunction();
                
        //         typename Nodes::VecLocD x;

        //         double p_int = 0.;
       
        //         if (fabs(nodesPFine_[connecPFine(i)] -> getPressure()) > 0){
        //             p_int = nodesPFine_[connecPFine(i)] -> getPressure();
        //         }else{
        //             p_int = nodesPFine_[connecPFine(i)] -> getPressure() + p;
        //         };

        //         nodesPFine_[connecPFine(i)] -> setPressureArlequin(p_int);
        //     };
        // };
        
        // for (int i=0; i<numNodesPCoarse; i++){
        //     nodesPCoarse_[i] -> 
        //         setPressureArlequin(nodesPCoarse_[i] -> getPressure());
        // };
        
        
        if (rank == 0) {
            //Printing results
            printVelocity(iTimeStep);
        };

    };
        
    return 0;

};



#endif
