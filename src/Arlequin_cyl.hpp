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

/// Mounts the overlapping mesh problem for solving the incompressible flow problem

template<int DIM>
class Arlequin{
public:
    /// Defines the class Fluid locally
    typedef Fluid<DIM>                     FluidMesh;

    /// Defines the class Element locally
    typedef typename FluidMesh::Elements   Elements;

    /// Defines the class Node locally
    typedef typename FluidMesh::Node       Nodes;

    /// Defines the class Boundary locally
    typedef typename FluidMesh::Boundaries Boundary;

    /// Defines the class SpecialQuad locally
    typedef typename Elements::SpecialQuad Quadrature;

    /// Defines the class Glue locally
    typedef Glue<DIM>                      GlueZone;

    typedef FluidParameters<DIM> Parameters;

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


private:
    int numElemCoarse;
    int numElemFine;
    int numBoundElemCoarse;
    int numBoundElemFine;
    int numElemGlueZoneFine;
    int numElemGlueZoneCoarse;
    int numNodesCoarse;
    int numNodesFine;
    int numNodesGlueZoneFine;
    int numNodesGlueZoneCoarse;
    int numTimeSteps;
    double dTime;
    int rank;
    int iTimeStep;

    Parameters *parametersCoarse, *parametersFine;

    std::pair<idx_t*,idx_t*> domDecompCoarse;//Coarse Model Domain Decomposition
    std::pair<idx_t*,idx_t*> domDecompFine;  //Fine Model Domain Decomposition

    Quadrature quad;

    double pi = M_PI;

    double alpha_f;
    double alpha_m;
    double gamma;

public:
    /// Sets the coarse and mesh models. It is considered that the fine model
    /// is completely immersed on the coarse model.
    /// @param Fluid coarse model @param Fluid fine model
    void setFluidModels(FluidMesh& coarse, FluidMesh& fine);

    /// Mounts and solve the incompressible flow problem with overlapping meshes
    /// using the Arlequin method whit the gluing zone defined in the fine model
    /// @param int maximum number of iterations of the Newton-Raphson's process
    /// @param double tolerance of the Newton-Raphson's process
    /// @param int type of problem to be solved: 1 - Stokes; 2 - Navier-Stokes.
    /// @param int 0 - Steady problem; 1 - Transient problem.
    int solveArlequinProblem(int iterNumber, double tolerance,
                             int problem_type, int time_dependency);

    /// Mounts and solve the incompressible flow problem with overlapping meshes
    /// using the Arlequin method with the gluing zone defined in the fine model
    /// and the fine model can be moved arbitrarily in an ALE description 
    /// framework
    /// @param int maximum number of iterations of the Newton-Raphson's process
    /// @param double tolerance of the Newton-Raphson's process
    /// @param int type of problem to be solved: 1 - Stokes; 2 - Navier-Stokes.
    /// @param int 0 - Steady problem; 1 - Transient problem.
    int solveArlequinProblemMoving(int iterNumber, double tolerance,
                                   int problem_type, int time_dependency);

    /// Mounts and solve the incompressible flow problem with overlapping meshes
    /// using the Arlequin method with the gluing zone defined in the fine model
    /// and the fine model can be moved arbitrarily in an ALE description 
    /// framework for the FSI interaction problem
    /// @param int maximum number of iterations of the Newton-Raphson's process
    /// @param double tolerance of the Newton-Raphson's process
    /// @param int type of problem to be solved: 1 - Stokes; 2 - Navier-Stokes.
    int solveFSIArlequin(int iterNumber, double tolerance,
                         int problem_type);

    /// Mounts and solve the incompressible flow problem with overlapping meshes
    /// using the Arlequin method whit the gluing zone defined in the coarse
    ///  model
    /// @param int maximum number of iterations of the Newton-Raphson's process
    /// @param double tolerance of the Newton-Raphson's process
    /// @param int type of problem to be solved: 1 - Stokes; 2 - Navier-Stokes.
    /// @param int 0 - Steady problem; 1 - Transient problem.
    int solveArlequinProblemCoarse(int iterNumber, double tolerance,
                                   int problem_type, int time_dependency);

    /// Compute and store the element boxes for improving the correspondence 
    /// searching process
    void setElementBoxes();

    /// Defines the gluing (or coupling) zone
    void setCouplingZone();

    /// Compute and store the signaled distance function
    void setSignaledDistance();

    /// Sets the energy weight function 
    /// @param double reference value for computing the energy weight function
    void setWeightFunction(double val);

    /// Compute the energy weight function for the coarse model
    /// @param double reference value (radius from the reference value)
    /// @param double epsilon for the Arlequin method to be relevant
    double weightFunctionCoarseValue(double r, double epsilon);

    /// Compute the energy weight function for the fine model
    /// @param double reference value (radius from the reference value)
    /// @param double epsilon for the Arlequin method to be relevant
    double weightFunctionFineValue(double r, double epsilon);

    /// Sets the nodal correspondence of the fine to the coarse models
    void setNodalCorrespondenceFine();

    /// Sets the nodal correspondence of the fine to the coarse models
    void setNodalCorrespondenceCoarse();

    /// Compute and print drag and lift coefficients
    void dragAndLiftCoefficients(std::ofstream& dragLift);

    /// Searchs the point correspondence in a fluid model
    /// @param VecLocD point
    /// @param vector<Nodes> vector of fluid model nodes
    /// @param vector<Elements> vector of fluid model elements
    /// @param int number of elements of the fluid model
    std::pair<int,ublas::bounded_vector<double,2> > 
                      searchNodeCorrespondence(typename Nodes::VecLocD x,
                                               std::vector<Nodes *> nodes,
                                               std::vector<Elements *> elements,
                                               int numElem);

    /// Print the results for Paraview post-processing
    /// @param int time step
    void printResults(int step);

    void initialAcceleration();

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
    ublas::bounded_vector<double,2> xk, Xk;

    //Compute element boxes for coarse model
    //Only function for straight elements
    for (int jel = 0; jel < numElemCoarse; jel++){
        connec = elementsCoarse_[jel] -> getConnectivity();
        x1 = nodesCoarse_[connec(0)] -> getCoordinates();
        x2 = nodesCoarse_[connec(1)] -> getCoordinates();
        x3 = nodesCoarse_[connec(2)] -> getCoordinates();      

        xk(0) = std::min(x1(0),std::min(x2(0), x3(0)));
        xk(1) = std::min(x1(1),std::min(x2(1), x3(1)));

        Xk(0) = std::max(x1(0),std::max(x2(0), x3(0)));
        Xk(1) = std::max(x1(1),std::max(x2(1), x3(1)));        
        
        elementsCoarse_[jel] -> setIntersectionParameters(xk, Xk);
    };

    //Compute element boxes for fine model
    //Only function for straight elements
    for (int jel = 0; jel < numElemFine; jel++){
        connec = elementsFine_[jel] -> getConnectivity();
        x1 = nodesFine_[connec(0)] -> getCoordinates();
        x2 = nodesFine_[connec(1)] -> getCoordinates();
        x3 = nodesFine_[connec(2)] -> getCoordinates();      

        xk(0) = std::min(x1(0),std::min(x2(0), x3(0)));
        xk(1) = std::min(x1(1),std::min(x2(1), x3(1)));

        Xk(0) = std::max(x1(0),std::max(x2(0), x3(0)));
        Xk(1) = std::max(x1(1),std::max(x2(1), x3(1)));        
        
        elementsFine_[jel] -> setIntersectionParameters(xk, Xk);
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
    //int numberIntPoints = elementsFine_[0] -> getNumberOfIntegrationPoints();
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
                    setIntegrationPointCorrespondence(i,corresp.first,
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
                    setIntegrationPointCorrespondence(i,corresp.first,
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
    typename Nodes::VecLocD x, x1, x2, n, test;
    typename Boundary::BoundConnect bconnec;
    double dist;

    for (int inode = 0 ; inode < numNodesFine; inode++){ 
        nodesFine_[inode] -> clearInnerNormal();
        nodesFine_[inode] -> setDistFunction(0.0);
    };
    for (int inode = 0 ; inode < numNodesCoarse; inode++){ 
        nodesCoarse_[inode] -> setDistFunction(0.0);
    };
    //approximate normal calculation

    for (int i = 0; i < numBoundElemFine; i++){
        if (boundaryFine_[i]->getConstrain(0) == 2){

            connec = elementsFine_[boundaryFine_[i] -> getElement()] -> 
                getConnectivity();
                            
            if (elementsFine_[boundaryFine_[i] -> getElement()] ->
                getElemSideInBoundary() == 0){
                bconnec(0) = connec(1);
                bconnec(1) = connec(2);
                bconnec(2) = connec(4);
            };
            if (elementsFine_[boundaryFine_[i] -> getElement()] ->
                getElemSideInBoundary() == 1){
                bconnec(0) = connec(2);
                bconnec(1) = connec(0);
                bconnec(2) = connec(5);
            };
            if (elementsFine_[boundaryFine_[i] -> getElement()] ->
                getElemSideInBoundary() == 2){
                bconnec(0) = connec(0);
                bconnec(1) = connec(1);
                bconnec(2) = connec(3);
            };
            
            //bconnec = boundaryFine_[i]->getBoundaryConnectivity();
            
            //first segment
            int no1 = bconnec(0);
            int no2 = bconnec(2);
            x1 = nodesFine_[no1] -> getCoordinates();
            x2 = nodesFine_[no2] -> getCoordinates();

            double sLength = sqrt((x2(1) - x1(1)) * (x2(1) - x1(1)) +
                                  (x1(0) - x2(0)) * (x1(0) - x2(0)));
            n(0) = (x2(1) - x1(1)) / sLength;
            n(1) = (x1(0) - x2(0)) / sLength;

            nodesFine_[no1] -> setInnerNormal(n);
            nodesFine_[no2] -> setInnerNormal(n);

            //second segment
            no1 = bconnec(2);
            no2 = bconnec(1);
            x1 = nodesFine_[no1] -> getCoordinates();
            x2 = nodesFine_[no2] -> getCoordinates();

            sLength = sqrt((x2(1) - x1(1)) * (x2(1) - x1(1)) +
                           (x1(0) - x2(0)) * (x1(0) - x2(0)));

            n(0) = (x2(1) - x1(1)) / sLength;
            n(1) = (x1(0) - x2(0)) / sLength;

            nodesFine_[no1] -> setInnerNormal(n);
            nodesFine_[no2] -> setInnerNormal(n);
            
        };
    };


    //Fine mesh nodes
    for (int ino = 0; ino < numNodesFine; ino++){
        x = nodesFine_[ino] -> getCoordinates();
        dist=10000000000000000000000000000.;
        
        for (int i = 0; i < numBoundElemFine; i++){
            if (boundaryFine_[i] -> getConstrain(0) == 2){
                connec = elementsFine_[boundaryFine_[i] -> getElement()] -> 
                    getConnectivity();
                
                if (elementsFine_[boundaryFine_[i] -> getElement()] ->
                    getElemSideInBoundary() == 0){
                    bconnec(0) = connec(1);
                    bconnec(1) = connec(2);
                    bconnec(2) = connec(4);
                };
                if (elementsFine_[boundaryFine_[i] -> getElement()] ->
                    getElemSideInBoundary() == 1){
                    bconnec(0) = connec(2);
                    bconnec(1) = connec(0);
                    bconnec(2) = connec(5);
                };
                if (elementsFine_[boundaryFine_[i] -> getElement()] ->
                    getElemSideInBoundary() == 2){
                    bconnec(0) = connec(0);
                    bconnec(1) = connec(1);
                    bconnec(2) = connec(3);
                };
                //bconnec = boundaryFine_[i] -> getBoundaryConnectivity();

                //first segment
                int no1 = bconnec(0);
                int no2 = bconnec(2);
                // std::cout<<no1<<" nos "<<no2<<std::endl;
              
                x1 = nodesFine_[no1] -> getCoordinates();
                x2 = nodesFine_[no2] -> getCoordinates();
                
                double aux0 =  sqrt((x2(1) - x1(1)) * (x2(1) - x1(1)) +
                                    (x2(0) - x1(0)) * (x2(0) - x1(0)));
                double aux1 = ((x(0) - x1(0)) * (x2(0) - x1(0))+
                               (x(1) - x1(1)) * (x2(1) - x1(1))) / aux0;
                double dist2 =-((x2(1) - x1(1)) * x(0) - 
                                (x2(0) - x1(0)) * x(1) +
                                x2(0) * x1(1) - x2(1) * x1(0)) / aux0;
                
                if (aux1 > aux0){
                    dist2 = sqrt((x2(1) - x(1)) * (x2(1) - x(1)) +
                                 (x2(0) - x(0)) * (x2(0) - x(0)));
                    //find signal
                    //side normal vector
                    n = nodesFine_[no2] -> getInnerNormal();
                    
                    test(0) = x(0) - x2(0);
                    test(1) = x(1) - x2(1);
                    double signaltest = inner_prod(n,test);
                    double signal = -1.;
                    
                    if (signaltest <= -0.001)signal = 1.;
                    
                    dist2 *= signal;
                };

                if (aux1 < 0.){
                    dist2 = sqrt((x(1) - x1(1)) * (x(1) - x1(1)) +
                                 (x(0) - x1(0)) * (x(0) - x1(0)));
                    //find signal
                    //side normal vector
                    n = nodesFine_[no1] -> getInnerNormal();
                    
                    test(0) = x(0) - x1(0);
                    test(1) = x(1) - x1(1);
                    double signaltest = inner_prod(n,test);
                    double signal = -1.;
                    
                    if (signaltest <= -0.001) signal = 1.;
                    
                    dist2 *= signal;
                };
                
                if (fabs(dist2) < fabs(dist)) dist = dist2;
                
                //second segment
                no1 = bconnec(2);
                no2 = bconnec(1);
                x1 = nodesFine_[no1] -> getCoordinates();
                x2 = nodesFine_[no2] -> getCoordinates();
                
                aux0 = sqrt((x2(1) - x1(1)) * (x2(1) - x1(1)) +
                            (x2(0) - x1(0)) * (x2(0) - x1(0)));
                aux1 = ((x(0) - x1(0)) * (x2(0) - x1(0)) + 
                        (x(1) - x1(1)) * (x2(1) - x1(1))) / aux0;
                dist2 = -((x2(1) - x1(1)) * x(0) - (x2(0) - x1(0)) * x(1) +
                          x2(0) * x1(1) - x2(1) * x1(0)) / aux0;

                if (aux1 > aux0){
                    dist2 = sqrt((x2(1) - x(1)) * (x2(1) - x(1)) +
                                 (x2(0) - x(0)) * (x2(0) - x(0)));
                    n = nodesFine_[no2] -> getInnerNormal();
                    
                    test(0) = x(0) - x2(0);
                    test(1) = x(1) - x2(1);
                    double signaltest = inner_prod(n,test);
                    double signal = -1.;
                    
                    if (signaltest <= -0.001)signal = 1.;
                    
                    dist2 *= signal;
                                       
                };

                if (aux1 < 0.){
                    dist2 = sqrt((x(1) - x1(1)) * (x(1) - x1(1)) +
                                 (x(0) - x1(0)) * (x(0) - x1(0)));
                    //find signal
                    //side normal vector
                    n = nodesFine_[no1] -> getInnerNormal();
                    
                    test(0) = x(0) - x1(0);
                    test(1) = x(1) - x1(1);
                    double signaltest = inner_prod(n,test);
                    double signal = -1.;
                    
                    if (signaltest <= -0.001)signal = 1.;
                    
                    dist2 *= signal;
                    
                };
                if (fabs(dist2) < fabs(dist)) dist = dist2;
            }; //if bf is the glue boundary
        }; //
    
        if(dist < 0) dist = 0;
        nodesFine_[ino] -> setDistFunction(dist);
     };

    //Coarse mesh
     for (int ino = 0; ino < numNodesCoarse; ino++){
        x = nodesCoarse_[ino] -> getCoordinates();
        dist=10000000000000000000000000000.;
        
        for (int i = 0; i < numBoundElemFine; i++){
            if (boundaryFine_[i] -> getConstrain(0) == 2){

                connec = elementsFine_[boundaryFine_[i] -> getElement()] -> 
                    getConnectivity();
                
                if (elementsFine_[boundaryFine_[i] -> getElement()] ->
                    getElemSideInBoundary() == 0){
                    bconnec(0) = connec(1);
                    bconnec(1) = connec(2);
                    bconnec(2) = connec(4);
                };
                if (elementsFine_[boundaryFine_[i] -> getElement()] ->
                    getElemSideInBoundary() == 1){
                    bconnec(0) = connec(2);
                    bconnec(1) = connec(0);
                    bconnec(2) = connec(5);
                };
                if (elementsFine_[boundaryFine_[i] -> getElement()] ->
                    getElemSideInBoundary() == 2){
                    bconnec(0) = connec(0);
                    bconnec(1) = connec(1);
                    bconnec(2) = connec(3);
                };
                //bconnec = boundaryFine_[i] -> getBoundaryConnectivity();

                //first segment
                int no1 = bconnec(0);
                int no2 = bconnec(2);
                // std::cout<<no1<<" nos "<<no2<<std::endl;
              
                x1 = nodesFine_[no1] -> getCoordinates();
                x2 = nodesFine_[no2] -> getCoordinates();
                
                double aux0 =  sqrt((x2(1) - x1(1)) * (x2(1) - x1(1)) +
                                    (x2(0) - x1(0)) * (x2(0) - x1(0)));
                double aux1 = ((x(0) - x1(0)) * (x2(0) - x1(0))+
                               (x(1) - x1(1)) * (x2(1) - x1(1))) / aux0;
                double dist2 =-((x2(1) - x1(1)) * x(0) - 
                                (x2(0) - x1(0)) * x(1) +
                                x2(0) * x1(1) - x2(1) * x1(0)) / aux0;
                
                if (aux1 > aux0){
                    dist2 = sqrt((x2(1) - x(1)) * (x2(1) - x(1)) +
                                 (x2(0) - x(0)) * (x2(0) - x(0)));
                    //find signal
                    //side normal vector
                    n = nodesFine_[no2] -> getInnerNormal();
                    
                    test(0) = x(0) - x2(0);
                    test(1) = x(1) - x2(1);
                    double signaltest = inner_prod(n,test);
                    double signal = -1.;
                    
                    if (signaltest <= -0.001)signal = 1.;
                    
                    dist2 *= signal;
                };

                if (aux1 < 0.){
                    dist2 = sqrt((x(1) - x1(1)) * (x(1) - x1(1)) +
                                 (x(0) - x1(0)) * (x(0) - x1(0)));
                    //find signal
                    //side normal vector
                    n = nodesFine_[no1] -> getInnerNormal();
                    
                    test(0) = x(0) - x1(0);
                    test(1) = x(1) - x1(1);
                    double signaltest = inner_prod(n,test);
                    double signal = -1.;
                    
                    if (signaltest <= -0.001) signal = 1.;
                    
                    dist2 *= signal;
                };
                
                if (fabs(dist2) < fabs(dist)) dist = dist2;
                
                //second segment
                no1 = bconnec(2);
                no2 = bconnec(1);
                x1 = nodesFine_[no1] -> getCoordinates();
                x2 = nodesFine_[no2] -> getCoordinates();
                
                aux0 = sqrt((x2(1) - x1(1)) * (x2(1) - x1(1)) +
                            (x2(0) - x1(0)) * (x2(0) - x1(0)));
                aux1 = ((x(0) - x1(0)) * (x2(0) - x1(0)) + 
                        (x(1) - x1(1)) * (x2(1) - x1(1))) / aux0;
                dist2 = -((x2(1) - x1(1)) * x(0) - (x2(0) - x1(0)) * x(1) +
                          x2(0) * x1(1) - x2(1) * x1(0)) / aux0;

                if (aux1 > aux0){
                    dist2 = sqrt((x2(1) - x(1)) * (x2(1) - x(1)) +
                                 (x2(0) - x(0)) * (x2(0) - x(0)));
                    n = nodesFine_[no2] -> getInnerNormal();
                    
                    test(0) = x(0) - x2(0);
                    test(1) = x(1) - x2(1);
                    double signaltest = inner_prod(n,test);
                    double signal = -1.;
                    
                    if (signaltest <= -0.001)signal = 1.;
                    
                    dist2 *= signal;
                                       
                };

                if (aux1 < 0.){
                    dist2 = sqrt((x(1) - x1(1)) * (x(1) - x1(1)) +  
                                 (x(0) - x1(0)) * (x(0) - x1(0)));
                    //find signal
                    //side normal vector
                    n = nodesFine_[no1] -> getInnerNormal();
                    
                    test(0) = x(0) - x1(0);
                    test(1) = x(1) - x1(1);
                    double signaltest = inner_prod(n,test);
                    double signal = -1.;
                    
                    if (signaltest <= -0.001)signal = 1.;
                    
                    dist2 *= signal;
                    
                };
                if (fabs(dist2) < fabs(dist)) dist = dist2;
            }; //if bf is the glue boundary
        }; //
    
        if (fabs(nodesCoarse_[ino] -> getDistFunction()) < 1.e-2){
            nodesCoarse_[ino] -> setDistFunction(dist); 
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
    int nodesCZ2[numNodesCoarse];

    double lim1 = 0.06251;
    double lim2 = 0.93749;
    double tick = 0.01249;

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
            double dist = nodesFine_[connec(ino)] -> getDistFunction();
            //  std::cout << "DIST " << dist << std::endl;
            if (dist <= fineModel.glueZoneThickness + 0.001){
                //if ((x(1) > 0.001) && (x(1) < 3.99)){
                  flag += 1;
                //    break;
              //};
            };
        };

        if (flag == 6) {
            elementsGlueZoneFine_.push_back(jel);
            elementsFine_[jel] -> setGlueZone();

            GlueZone *el = new GlueZone(index++,jel);
            glueZoneFine_.push_back(el);
            
            // for (int i=0; 
            //      i < elementsFine_[jel] -> getNumberOfIntegrationPoints(); i++){
                
            //     x = elementsFine_[jel] -> getIntegPointCoordinatesValue(i);  
                
            //     //double dist = nodesFine_[connec(ino)] -> getDistFunction();
            //     double dist = 0;
            //     if (dist <= 1.01){
            //         elementsFine_[jel] -> setIntegPointInGlueZone(i);    
            //     };

            //     // if ((x(0) < lim1) || (x(0) > lim2) || 
            //     //     (x(1) < lim1) || (x(1) > lim2)){
            //     // if ((x(0) < 4.) || (x(0) > 12.5) || 
            //     //     (x(1) < 4.) || (x(1) > 8.)){
            //     //     elementsFine_[jel] -> setIntegPointInGlueZone(i);
            //     // };

            // };
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

    if (rank == 0) std::cout << "GLUE ZONE - Number of Nodes = " 
                             << numNodesGlueZoneFine 
                             << " - Number of Elements = " 
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

    // Glue Zone in coarse mesh
    for (int i = 0; i < numNodesCoarse; i++) nodesCZ2[i] = 0;    

    elementsGlueZoneCoarse_.reserve(numElemCoarse / 3);
    nodesGlueZoneCoarse_.reserve(numNodesCoarse / 3);

    //Defines a criterion to select the elements that are in the glue zone
    for (int jel = 0; jel < numElemCoarse; jel++){
        
        connec = elementsCoarse_[jel] -> getConnectivity();
        flag = 0;

        for (int ino = 0; ino < 6; ino++){
            x = nodesCoarse_[connec(ino)] -> getCoordinates();

            if ((x(0) < lim1) || (x(0) > lim2) || 
                (x(1) < lim1) || (x(1) > lim2)){

            }else{
                if ((x(0) > lim1 + tick) && (x(0) < lim2 - tick) &&
                    (x(1) > lim1 + tick) && (x(1) < lim2 - tick)){

                }else{
                    flag = 1;
                    break;
                };
            };
        };
        if (flag > 0) {
            elementsGlueZoneCoarse_.push_back(jel);
            elementsCoarse_[jel] -> setGlueZone();

            GlueZone *el = new GlueZone(index++,jel);
            glueZoneCoarse_.push_back(el);
            
            for (int i=0; i < elementsCoarse_[jel] -> 
                     getNumberOfIntegrationPoints(); i++){
                
                x = elementsCoarse_[jel] -> getIntegPointCoordinatesValue(i);

                if ((x(0) < lim1) || (x(0) > lim2) || 
                    (x(1) < lim1) || (x(1) > lim2)){

                }else{
                    if ((x(0) > lim1 + tick) && (x(0) < lim2 - tick) &&
                        (x(1) > lim1 + tick) && (x(1) < lim2 - tick)){

                    }else{
                    elementsCoarse_[jel] -> setIntegPointInGlueZone(i);
                    };
                };
            };
        };        
    };

    //Defines which nodes are in the glue zone
    numElemGlueZoneCoarse = elementsGlueZoneCoarse_.size();
    for (int i = 0; i < numElemGlueZoneCoarse; i++){
        connec = elementsCoarse_[elementsGlueZoneCoarse_[i]] -> 
            getConnectivity();

        for (int ino = 0; ino < 6; ino++){
            nodesCZ2[connec(ino)] += 1;
        };
    };

    //Compute number of nodes in the glue zone
    numNodesGlueZoneCoarse = 0;
    for (int i = 0; i < numNodesCoarse; i++){
        if(nodesCZ2[i] > 0) {
            numNodesGlueZoneCoarse += 1;
            nodesGlueZoneCoarse_.push_back(i);
        };
    };

    for (int i = 0; i < numNodesGlueZoneCoarse; i++){
        typename Nodes::VecLocD x;
        x = nodesCoarse_[nodesGlueZoneCoarse_[i]] -> getCoordinates();
        
        Nodes *no = new Nodes(x,i);
        nodesLagrangeCoarse_.push_back(no);
    };

    for (int i = 0; i < numElemGlueZoneCoarse; i++){
        typename Elements::Connectivity connecAux;

        connec = elementsCoarse_[elementsGlueZoneCoarse_[i]] -> 
            getConnectivity();
        
        for (int ino = 0; ino < numNodesGlueZoneCoarse; ino++){
            if (nodesGlueZoneCoarse_[ino] == connec(0)) connecAux(0) = ino;
            if (nodesGlueZoneCoarse_[ino] == connec(1)) connecAux(1) = ino;
            if (nodesGlueZoneCoarse_[ino] == connec(2)) connecAux(2) = ino;
            if (nodesGlueZoneCoarse_[ino] == connec(3)) connecAux(3) = ino;
            if (nodesGlueZoneCoarse_[ino] == connec(4)) connecAux(4) = ino;
            if (nodesGlueZoneCoarse_[ino] == connec(5)) connecAux(5) = ino;
        };
        
        glueZoneCoarse_[i] -> setConnectivity(connecAux);
        glueZoneCoarse_[i] -> setNodes(nodesLagrangeCoarse_);

    };
     
};

//------------------------------------------------------------------------------
//----------------------COMPUTES THE WEIGHT FUNCTION VALUE----------------------
//------------------------------------------------------------------------------
template<>
double Arlequin<2>::weightFunctionFineValue(double r, double epsilon){

    double wFuncValue;

    if (r <= 4.5){
        wFuncValue = 1. - epsilon;
    } else {
        if (r > 5.55){
            wFuncValue = 0.;
        } else {
            //wFuncValue = 1. - epsilon;

            //Linear weight function
            // 0.5 < r < 3.5
            //wFuncValue = - (1. - epsilon) * (r - 3.5) / 3.;
            // 2.5 < r < 3.5
            wFuncValue = - (1. - epsilon) * (r - 3.5);
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

    if (r <= 4.5){
        wFuncValue = epsilon;
    } else {
        if (r >= 5.5){
            wFuncValue = 1.;
        } else {
            //wFuncValue = epsilon;

            //Linear weight function
            // 0.5 < r < 3.5
            //wFuncValue = (1. - epsilon) * (r - 3.5) / 3. + 1.;
            // 2.5 < r < 3.5
            wFuncValue = (1. - epsilon) * (r - 3.5) + 1.;
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
    double wFuncValue;

    double epsilon = coarseModel.arlequinEpsilon;
    double lambda = fineModel.glueZoneThickness*1.01;
 
    for (int i = 0; i < numNodesCoarse; i++){
        
        double r = nodesCoarse_[i] -> getDistFunction();
            
        if (r < 0){
            wFuncValue = 1.;
            //if (fabs(r) < 1.e-5) wFuncValue = 0.5;
        } else {
            if (r >= lambda){
                wFuncValue = epsilon;
            } else {
                wFuncValue = 1. - (1. - epsilon) / lambda * r;
                //wFuncValue = 0.5;
                //wFuncValue = epsilon;
                // wFuncValue = 1+3.*(epsilon-1.)/(lambda*lambda) * r * r
                //     -2.*(epsilon-1.)/(lambda*lambda*lambda) * r * r * r;
                
                if (wFuncValue < epsilon) wFuncValue = epsilon;
            };
        };  
      
        //wFuncValue = weightFunctionCoarseValue(r,epsilon);
        
        nodesCoarse_[i] -> setWeightFunction(wFuncValue);
    };         

    for (int jel = 0; jel < numElemCoarse; jel++){
        elementsCoarse_[jel] -> setIntegPointWeightFunction();        
    };


    for (int i=0; i<numNodesFine; i++){

        double r = nodesFine_[i] -> getDistFunction();
        
        if (r >= lambda){
            wFuncValue = 1. - epsilon;
        } else {
            wFuncValue = (1. - epsilon) / lambda * r;
            //wFuncValue = 0.5;
            //wFuncValue = 1. - epsilon;
            // wFuncValue = -3.*(epsilon-1.)/(lambda*lambda) * r * r
            //         +2.*(epsilon-1.)/(lambda*lambda*lambda) * r * r * r;

            if (wFuncValue > (1. - epsilon)) wFuncValue = 1. - epsilon;

        };  
        
        // wFuncValue = weightFunctionFineValue(r,epsilon);

        nodesFine_[i] -> setWeightFunction(wFuncValue);
    };
    for (int jel = 0; jel < numElemFine; jel++){
        elementsFine_[jel] -> setIntegPointWeightFunction();        
    };

    

     
    return;
};

//------------------------------------------------------------------------------
//----------------------------PRINT VELOCITY RESULTS----------------------------
//------------------------------------------------------------------------------
template<>
void Arlequin<2>::printResults(int step) {

    //PRINT COARSE MODEL RESULTS
    
    std::string result;
    std::ostringstream convert;

    convert << step+100000;
    result = convert.str();
    std::string s = "saidaVelCoarse"+result+".vtu";
    
    std::fstream output_v(s.c_str(), std::ios_base::out);

    output_v << "<?xml version=\"1.0\"?>" << std::endl
             << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl
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
        std::string b;
        std::ostringstream a;
        a << x(0) << " " << x(1) << " " << 0.0 << " ";
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

    if (coarseModel.printVelocity){
        output_v <<"      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                 << "Name=\"Velocity\" format=\"ascii\">" << std::endl;
        for (int i=0; i<numNodesCoarse; i++){        
            output_v << nodesCoarse_[i] -> getVelocity(0) << " "              
                     << nodesCoarse_[i] -> getVelocity(1) << " " 
                     << 0. << std::endl;
        };
        output_v << "      </DataArray> " << std::endl;
        output_v <<"      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                 << "Name=\"Acceleration\" format=\"ascii\">" << std::endl;
        for (int i=0; i<numNodesCoarse; i++){        
            output_v << nodesCoarse_[i] -> getAcceleration(0) << " "              
                     << nodesCoarse_[i] -> getAcceleration(1) << " " 
                     << 0. << std::endl;
        };
        output_v << "      </DataArray> " << std::endl;
    };

    if (coarseModel.printRealVelocity){
        output_v <<"      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                 << "Name=\"Real Velocity\" format=\"ascii\">" << std::endl;
        for (int i=0; i<numNodesCoarse; i++){        
            output_v << nodesCoarse_[i] -> getVelocityArlequin(0) << " " 
                     << nodesCoarse_[i] -> getVelocityArlequin(1) << " " 
                     << 0. << std::endl;
        };
        output_v << "      </DataArray> " << std::endl;
    }

    if (coarseModel.printLagrangeMultipliers){
        output_v <<"      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                 << "Name=\"Lagrange Multipliers\" format=\"ascii\">"
                 << std::endl;
        for (int i=0; i<numNodesCoarse; i++){        
            output_v << nodesCoarse_[i] -> getLagrangeMultiplier(0) << " "
                     << nodesCoarse_[i] -> getLagrangeMultiplier(1) << " " 
                     << 0. << std::endl;
        };
        output_v << "      </DataArray> " << std::endl;
    }

    if (coarseModel.printElementCorrespondence){
        output_v <<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                 << "Name=\"Element\" format=\"ascii\">" << std::endl;
        for (int i=0; i<numNodesCoarse; i++){
            output_v << nodesCoarse_[i] -> getNodalElemCorrespondence()
                     << std::endl;
        };
        output_v << "      </DataArray> " << std::endl;
    }
    
    if (coarseModel.printDistFunction){
        output_v <<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                 << "Name=\"Dist Function\" format=\"ascii\">" << std::endl;
        for (int i=0; i<numNodesCoarse; i++){        
            output_v << nodesCoarse_[i] -> getDistFunction() << std::endl;
        };
        output_v << "      </DataArray> " << std::endl;
    }

    if (coarseModel.printEnergyWeightFunction){
        output_v<< "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                << "Name=\"Weight Function\" format=\"ascii\">" << std::endl;
        for (int i=0; i<numNodesCoarse; i++){
            output_v << nodesCoarse_[i] -> getWeightFunction() << std::endl;
        };
        output_v << "      </DataArray> " << std::endl;
    }

    if (coarseModel.printPressure){
        output_v <<"      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                 << "Name=\"Pressure\" format=\"ascii\">" << std::endl;
        for (int i=0; i<numNodesCoarse; i++){
            output_v << 0. << " " << 0. << " " 
                     << nodesCoarse_[i] -> getPressure() << std::endl;
        };
        output_v << "      </DataArray> " << std::endl;
    }

    if (coarseModel.printRealPressure){
        output_v <<"      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                 << "Name=\"Real Pressure\" format=\"ascii\">" << std::endl;
        for (int i=0; i<numNodesCoarse; i++){
            output_v << 0. << " " << 0. << " " 
                     << nodesCoarse_[i] -> getPressureArlequin() << std::endl;
        };
        output_v << "      </DataArray> " << std::endl;
    };

    if (coarseModel.printVorticity){
        output_v <<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                 << "Name=\"Vorticity\" format=\"ascii\">" << std::endl;
        for (int i=0; i<numNodesCoarse; i++){
            output_v << nodesCoarse_[i] -> getVorticity() << std::endl;
        };
        output_v << "      </DataArray> " << std::endl;
    };


    output_v << "    </PointData>" << std::endl; 

    //WRITE ELEMENT RESULTS
    output_v << "    <CellData>" << std::endl;
    
    if (coarseModel.printProcess){
        output_v <<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                 << "Name=\"Process\" format=\"ascii\">" << std::endl;
        for (int i=0; i<numElemCoarse; i++){
            output_v << domDecompCoarse.first[i] << std::endl;
        };
        output_v << "      </DataArray> " << std::endl;
    }

    if (coarseModel.printEnergyWeightFunction){
        output_v <<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                 << "Name=\"Weight Function\" format=\"ascii\">" << std::endl;
        for (int i=0; i<numElemCoarse; i++){
            output_v << elementsCoarse_[i] -> getIntegPointWeightFunction(0)
                     << std::endl;
        };
        output_v << "      </DataArray> " << std::endl;
    }

    output_v <<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
             << "Name=\"Lines\" format=\"ascii\">" << std::endl;
    for (int i=0; i<numElemCoarse; i++){
        int res = 0;
        for (int j=0; j<numBoundElemCoarse; j++){
           if (boundaryCoarse_[j] -> getElement() == i) res = boundaryCoarse_[j] -> getBoundaryGroup();
        }
        output_v << res << std::endl;
    };
    output_v << "      </DataArray> " << std::endl;


    int cont=0;
    
    if (coarseModel.printGlueZone){
        output_v <<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                 << "Name=\"Glue Zone\" format=\"ascii\">" << std::endl;
        cont = 0;
        for (int i=0; i<numElemCoarse; i++){
            if (elementsGlueZoneCoarse_[cont] == i){
                output_v << 1.0 << std::endl;
                cont++;
            }else{
                output_v << 0.0 << std::endl;
            };
        };
        output_v << "      </DataArray> " << std::endl;
    }

    if (coarseModel.printJacobian){
        output_v <<"      <DataArray type=\"Float32\" NumberOfComponents=\"1\" "
                 << "Name=\"Jacobian\" format=\"ascii\">" << std::endl;
        for (int i=0; i<numElemCoarse; i++){
            output_v << elementsCoarse_[i] -> getJacobian() << std::endl;
        };
        output_v << "      </DataArray> " << std::endl;
    }


    // std::cout << lala << std::endl;

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
        output_vf << x(0) << " " << x(1) << " " << 0.001 << std::endl;        
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

    if (fineModel.printVelocity){
        output_vf<<"      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                  << "Name=\"Velocity\" format=\"ascii\">" << std::endl;
        for (int i=0; i<numNodesFine; i++){
            output_vf << nodesFine_[i] -> getVelocity(0) << " "              
                      << nodesFine_[i] -> getVelocity(1) << " " 
                      << 0. << std::endl;
        };
        output_vf << "      </DataArray> " << std::endl;
        output_vf<<"      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                  << "Name=\"Acceleration\" format=\"ascii\">" << std::endl;
        for (int i=0; i<numNodesFine; i++){
            output_vf << nodesFine_[i] -> getAcceleration(0) << " "              
                      << nodesFine_[i] -> getAcceleration(1) << " " 
                      << 0. << std::endl;
        };
        output_vf << "      </DataArray> " << std::endl;
    }

    if (fineModel.printInnerNormal){
        output_vf<<"      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                 << "Name=\"Inner Normal\" format=\"ascii\">" << std::endl;
        typename Nodes::VecLocD n;
        for (int i=0; i<numNodesFine; i++){
            n = nodesFine_[i] -> getInnerNormal();
            output_vf << n(0) << " " << n(1) << " " << 0. << std::endl;
        };
        output_vf << "      </DataArray> " << std::endl;
    }

    if (fineModel.printRealVelocity){
        output_vf<<"      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                 << "Name=\"Real Velocity\" format=\"ascii\">" << std::endl;
        for (int i=0; i<numNodesFine; i++){
            output_vf << nodesFine_[i] -> getVelocityArlequin(0) << " "
                      << nodesFine_[i] -> getVelocityArlequin(1) << " " 
                      << 0. << std::endl;
        };
        output_vf << "      </DataArray> " << std::endl;
    }

    if (fineModel.printMeshVelocity){
        output_vf<<"      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                 << "Name=\"Mesh Velocity\" format=\"ascii\">" << std::endl;
        for (int i=0; i<numNodesFine; i++){
            output_vf << nodesFine_[i] -> getMeshVelocity(0) << " " 
                      << nodesFine_[i] -> getMeshVelocity(1) << " " 
                      << 0. << std::endl;
        };
        output_vf << "      </DataArray> " << std::endl;
    };

    if (fineModel.printLagrangeMultipliers){
        output_vf<<"      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                 << "Name=\"Lagrange Multipliers\" format=\"ascii\">" << std::endl;
        for (int i=0; i<numNodesFine; i++){
            output_vf << nodesFine_[i] -> getLagrangeMultiplier(0) << " " 
                      << nodesFine_[i] -> getLagrangeMultiplier(1) << " " 
                      << 0. << std::endl;
        };
        output_vf << "      </DataArray> " << std::endl;
    }

    if (fineModel.printDistFunction){
        output_vf<<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                 << "Name=\"Dist Function\" format=\"ascii\">" << std::endl;
        for (int i=0; i<numNodesFine; i++){
            output_vf << nodesFine_[i] -> getDistFunction() << std::endl;
        };
        output_vf << "      </DataArray> " << std::endl;
    }

    if (fineModel.printElementCorrespondence){
        output_vf<<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                 << "Name=\"Element\" format=\"ascii\">" << std::endl;
        for (int i=0; i<numNodesFine; i++){
            output_vf << nodesFine_[i] -> getNodalElemCorrespondence() 
                      << std::endl;
        };
        output_vf << "      </DataArray> " << std::endl;
    }

    if (fineModel.printEnergyWeightFunction){
        output_vf<<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                 << "Name=\"Weight Function\" format=\"ascii\">" << std::endl;
        for (int i=0; i<numNodesFine; i++){
            output_vf << nodesFine_[i] -> getWeightFunction() << std::endl;
        };
        output_vf << "      </DataArray> " << std::endl;
    }
    
    if (fineModel.printMeshDisplacement){
        output_vf<<"      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                 << "Name=\"Mesh Displacement\" format=\"ascii\">" << std::endl;
        for (int i=0; i<numNodesFine; i++){
            typename Nodes::VecLocD x, xi;
            x = nodesFine_[i] -> getCoordinates();
            xi = nodesFine_[i] -> getInitialCoordinates();
            output_vf << x(0) - xi(0) << " " << x(1) - xi(1) << " " << 0. 
                      << std::endl;
        };
        output_vf << "      </DataArray> " << std::endl;
    };

    if (fineModel.printPressure){
        output_vf<<"      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                 << "Name=\"Pressure\" format=\"ascii\">" << std::endl;
        for (int i=0; i<numNodesFine; i++){
            output_vf << 0. << " " << 0. << " " 
                      << nodesFine_[i] -> getPressure() << std::endl;
        };
        output_vf << "      </DataArray> " << std::endl;
    };

    if (fineModel.printRealPressure){
        output_vf<<"      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                 << "Name=\"Real Pressure\" format=\"ascii\">" << std::endl;
        for (int i=0; i<numNodesFine; i++){
            output_vf << 0. << " " << 0. << " " 
                      << nodesFine_[i] -> getPressureArlequin() << std::endl;
        };
        output_vf << "      </DataArray> " << std::endl;
    }

    if (fineModel.printVorticity){
        output_vf<<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                 << "Name=\"Vorticity\" format=\"ascii\">" << std::endl;
        for (int i=0; i<numNodesFine; i++){
            output_vf << nodesFine_[i] -> getVorticity() << std::endl;
        };
        output_vf << "      </DataArray> " << std::endl;
    }

    output_vf<<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
             << "Name=\"Element\" format=\"ascii\">" << std::endl;
    for (int i=0; i<numNodesFine; i++){
        output_vf << nodesFine_[i] -> getNodalElemCorrespondence() 
                  << std::endl;
    };
    output_vf << "      </DataArray> " << std::endl;


    output_vf << "    </PointData>" << std::endl; 

    //WRITE ELEMENT RESULTS
    output_vf << "    <CellData>" << std::endl;
    
    if (fineModel.printProcess){
        output_vf<<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                 << "Name=\"Process\" format=\"ascii\">" << std::endl;
        for (int i=0; i<numElemFine; i++){
            output_vf << domDecompFine.first[i] << std::endl;
        };
        output_vf << "      </DataArray> " << std::endl;
    }

    output_vf <<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
             << "Name=\"Lines\" format=\"ascii\">" << std::endl;
    for (int i=0; i<numElemFine; i++){
        int res = 0;
        for (int j=0; j<numBoundElemFine; j++){
           if (boundaryFine_[j] -> getElement() == i) res = boundaryFine_[j] -> getBoundaryGroup();
        }
        output_vf << res << std::endl;
    };
    output_vf << "      </DataArray> " << std::endl;


    if (fineModel.printGlueZone){
        output_vf<<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
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
    }


    if (fineModel.printEnergyWeightFunction){
        output_vf<<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                 << "Name=\"Weight Function\" format=\"ascii\">" << std::endl;
        for (int i=0; i<numElemFine; i++){
            output_vf << elementsFine_[i] -> getIntegPointWeightFunction(0)
                      << std::endl;
        };
        output_vf << "      </DataArray> " << std::endl;
    };

    if (fineModel.printJacobian){
        output_vf<<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                 << "Name=\"Jacobian\" format=\"ascii\">" << std::endl;
        for (int i=0; i<numElemFine; i++){
            output_vf << elementsFine_[i] -> getJacobian() << std::endl;
        };
        output_vf << "      </DataArray> " << std::endl;
    };

    output_vf << "    </CellData>" << std::endl; 

    //FINALIZE OUTPUT FILE
    output_vf << "  </Piece>" << std::endl
           << "  </UnstructuredGrid>" << std::endl
           << "</VTKFile>" << std::endl;


    
    
};


//------------------------------------------------------------------------------
//------------SETS COARSE/FINE MESHES AND GETS ITS BASIC INFORMATIONS-----------
//------------------------------------------------------------------------------
template<>
void Arlequin<2>::setFluidModels(FluidMesh& coarse, FluidMesh& fine){

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

    numBoundElemFine = boundaryFine_.size();
    numBoundElemCoarse = boundaryCoarse_.size();

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



    // Rotating fine mesh
    for (int i = 0; i < numNodesFine; ++i){
        typename Nodes::VecLocD x, xn;
        x = nodesFine_[i] -> getCoordinates();       
    
        double a = -10 * pi / 180;

        xn(0) = 0.5 + (x(0)-0.5) * cos(a) - (x(1)-0.0) * sin(a);
        xn(1) = 0.0 + (x(0)-0.5) * sin(a) + (x(1)-0.0) * cos(a);

        nodesFine_[i] -> setCoordinates(xn);
        nodesFine_[i] -> setPreviousCoordinates(0,xn(0));
        nodesFine_[i] -> setPreviousCoordinates(1,xn(1));
    }

    parametersFine = &fineModel.fluidParameters;
    parametersCoarse = &coarseModel.fluidParameters;







    setSignaledDistance();

    //Construct the glue zone based on some defined criterion
    setCouplingZone();

    //Computes the Weight function for all the finite elements
    setWeightFunction(16.); 


    if(rank == 0){
        std::cout << "---------------------ARLEQUIN DATA---------------------" 
                  << std::endl;
        std::cout << "Coarse Model: " << numNodesCoarse << " nodes, " << numElemCoarse << " elements." << std::endl;
        std::cout << "Fine Model: " << numNodesFine << " nodes, " << numElemFine << " elements." << std::endl;
        std::cout << "Lagrange Multipliers: " << numNodesGlueZoneFine << " nodes, " << numElemGlueZoneFine << " elements." << std::endl;
        std::cout << "Number of Degrees of Freedom: " << 3*numNodesCoarse + 3*numNodesFine + 2*numNodesGlueZoneFine << std::endl; 
    }

    //if(rank == 0) printResults(0);


    // nodesFine_[269] -> setWeightFunction(0.);
    // nodesFine_[80] -> setWeightFunction(0.);





};

//------------------------------------------------------------------------------
//----------------------COMPUTES DRAG AND LIFT COEFFICIENTS---------------------
//------------------------------------------------------------------------------
template<>
void Arlequin<2>::dragAndLiftCoefficients(std::ofstream& dragLift){



    double dragCoefficient = 0.;
    double liftCoefficient = 0.;
    double pressureDragCoefficient = 0.;
    double pressureLiftCoefficient = 0.;
    double frictionDragCoefficient = 0.;
    double frictionLiftCoefficient = 0.;
    double theta = 0.;
    double pitchingMomentCoefficient = 0.;
    double pMom = 0.;
    double per = 0.;

    double rhoInf = 1.0;
    double velocityInf[2];
    velocityInf[0] = 1;
    velocityInf[1] = 0.;

    for (int jel = 0; jel < numBoundElemFine; jel++){   
        

        
        double dForce = 0.;
        double lForce = 0.;
        double pDForce = 0.;
        double pLForce = 0.;
        double fDForce = 0.;
        double fLForce = 0.;



        for (int i=0; i<fineModel.numberOfLines; i++){
            //std::cout << "Bound group " << boundary_[jel] -> getBoundaryGroup() << std::endl;
            if (boundaryFine_[jel] -> getBoundaryGroup() == fineModel.dragAndLiftBoundary[i]){
                //std::cout << "line " << i << " " << jel << std::endl;
                int iel = boundaryFine_[jel] -> getElement();
                elementsFine_[iel] -> computeDragAndLiftForces();

                // double theta2 = elementsFine_[iel] -> computeSeparationAngle();
                // if (theta2>100) theta = theta2;

                //if ((theta > 0) && (rank == 0)) std::cout << "AQUI SEPARATION ANGLE " << iel << " " << theta << std::endl;

                pDForce = elementsFine_[iel] -> getPressureDragForce();
                pLForce = elementsFine_[iel] -> getPressureLiftForce();
                fDForce = elementsFine_[iel] -> getFrictionDragForce();
                fLForce = elementsFine_[iel] -> getFrictionLiftForce();
                dForce = elementsFine_[iel] -> getDragForce();
                lForce = elementsFine_[iel] -> getLiftForce();
                pMom += elementsFine_[iel] -> getPitchingMoment();
                per += elementsFine_[iel] -> getPerimeter();
            };
        };
        
        pressureDragCoefficient += pDForce ;/// 
            //(0.5 * rhoInf * velocityInf[0] * velocityInf[0]);
        pressureLiftCoefficient += pLForce ;/// 
            //(0.5 * rhoInf * velocityInf[0] * velocityInf[0]);
        
        frictionDragCoefficient += fDForce / 
            (0.5 * rhoInf * velocityInf[0] * velocityInf[0] * 1);
        frictionLiftCoefficient += fLForce / 
            (0.5 * rhoInf * velocityInf[0] * velocityInf[0] * 1);
        
        dragCoefficient += dForce / 
            (0.5 * rhoInf * velocityInf[0] * velocityInf[0] * 1);
        liftCoefficient += lForce / 
            (0.5 * rhoInf * velocityInf[0] * velocityInf[0] * 1);
        
    };
    
    pitchingMomentCoefficient = pMom / (rhoInf * velocityInf[0] * velocityInf[0] * per);

    if (rank == 0) {
        const int timeWidth = 13;
        const int numWidth = 13;
        dragLift << std::setprecision(5) << std::scientific;
        dragLift << std::left << std::setw(timeWidth) << iTimeStep * dTime;
        dragLift << std::setw(numWidth) << pressureDragCoefficient;
        dragLift << std::setw(numWidth) << pressureLiftCoefficient;
        dragLift << std::setw(numWidth) << frictionDragCoefficient;
        dragLift << std::setw(numWidth) << frictionLiftCoefficient;
        dragLift << std::setw(numWidth) << dragCoefficient;
        dragLift << std::setw(numWidth) << liftCoefficient;
        dragLift << std::setw(numWidth) << pitchingMomentCoefficient;
        //dragLift << std::setw(numWidth) << theta;
        dragLift << std::endl;
    }

 

}



//------------------------------------------------------------------------------
//----------------COMPUTE ARLEQUIN COUPLED NAVIER-STOKES PROBLEM----------------
//------------------------------------------------------------------------------
template<>
void Arlequin<2>::initialAcceleration(){

    Mat               A, F;
    Vec               b, u, All;
    PetscErrorCode    ierr;
    PetscInt          Istart, Iend, Ii, Ione, iterations;
    KSP               ksp; 
    PC                pc;
    VecScatter        ctx;
    PetscScalar       val;
    //MatNullSpace      nullsp;
    
    iTimeStep = 0.;

    //Computes the Nodal correspondence between fine nodes and coarse elements
    setNodalCorrespondenceFine();

    // Computes the system size
    int sysSize = 3 * numNodesCoarse + 3 * numNodesFine + 2 * numNodesGlueZoneFine;

    double integScheme = fineModel.integScheme;

    alpha_f = 1. / (1. + integScheme);
    alpha_m = 0.5 * (3. - integScheme) / (1. + integScheme);
    gamma = 0.5 + alpha_m - alpha_f;

    if (rank == 0) {std::cout << "--- INITIAL ACCELERATION ---" << std::endl;}
            
    boost::posix_time::ptime t1 =
        boost::posix_time::microsec_clock::local_time();
    
    // Preallocates the matrix
    ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                        sysSize, sysSize, 300, NULL, 900, NULL, &A); 
    
    //Create PETSc vectors
    ierr = VecCreate(PETSC_COMM_WORLD, &b); 
    ierr = VecSetSizes(b, PETSC_DECIDE, sysSize); 
    ierr = VecSetFromOptions(b);
    ierr = VecDuplicate(b, &u);
    ierr = VecDuplicate(b, &All);
                
    for (int i=0; i<sysSize; i++){
        double val = 1.e-15;
        ierr = MatSetValues(A,1,&i,1,&i,&val,ADD_VALUES);
    }
                                

    //------------------------------------------------------------------
    //------------------BEGIN OF LINEAR SYSTEM ASSEMBLY-----------------
    //------------------------------------------------------------------
    //Coarse mesh
    for (int jel = 0; jel < numElemCoarse; jel++){   
          
        if (domDecompCoarse.first[jel] == rank) {
              
            //Compute Element matrix
            elementsCoarse_[jel] -> getTransientNavierStokesInitial();
            
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
                        ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&Ajac(2*i  ,2*j  ),ADD_VALUES);
                    };
                    if (fabs(Ajac(2*i+1,2*j  )) >= 1.e-15){
                        int dof_i = 2 * connec(i) + 1;
                        int dof_j = 2 * connec(j);
                        ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&Ajac(2*i+1,2*j  ),ADD_VALUES);
                    };
                    if (fabs(Ajac(2*i  ,2*j+1)) >= 1.e-15){
                        int dof_i = 2 * connec(i);
                        int dof_j = 2 * connec(j) + 1;
                        ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&Ajac(2*i  ,2*j+1),ADD_VALUES);
                    };
                    if (fabs(Ajac(2*i+1,2*j+1)) >= 1.e-15){
                        int dof_i = 2 * connec(i) + 1;
                        int dof_j = 2 * connec(j) + 1;
                        ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&Ajac(2*i+1,2*j+1),ADD_VALUES);
                    };
                    //Matrix Q and Qt
                    if (fabs(Ajac(2*i  ,12+j)) >= 1.e-15){
                        int dof_i = 2 * connec(i);
                        int dof_j = 2 * numNodesCoarse + connec(j);
                        ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&Ajac(2*i  ,12+j),ADD_VALUES);
                    };                  
                    if (fabs(Ajac(12+j,2*i  )) >= 1.e-15){
                        int dof_i = 2 * connec(i);
                        int dof_j = 2 * numNodesCoarse + connec(j);
                        ierr = MatSetValues(A,1,&dof_j,1,&dof_i,&Ajac(12+j,2*i  ),ADD_VALUES);
                    };
                    if (fabs(Ajac(2*i+1,12+j)) >= 1.e-15){
                        int dof_i = 2 * connec(i) + 1;
                        int dof_j = 2 * numNodesCoarse + connec(j);
                        ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&Ajac(2*i+1,12+j),ADD_VALUES);
                    };
                    if (fabs(Ajac(12+j,2*i+1)) >= 1.e-15){
                        int dof_i = 2 * connec(i) + 1;
                        int dof_j = 2 * numNodesCoarse + connec(j);
                        ierr = MatSetValues(A,1,&dof_j,1,&dof_i,&Ajac(12+j,2*i+1),ADD_VALUES);
                    };
                    if (fabs(Ajac(12+i,12+j)) >= 1.e-15){
                        int dof_i = 2 * numNodesCoarse + connec(i);
                        int dof_j = 2 * numNodesCoarse + connec(j);
                        ierr = MatSetValues(A,1,&dof_j,1,&dof_i,&Ajac(12+i,12+j),ADD_VALUES);
                    };
                };                  
                //Rhs vector
                if (fabs(Rhs(2*i  )) >= 1.e-15){
                    int dof_i = 2 * connec(i);
                    ierr = VecSetValues(b,1,&dof_i,&Rhs(2*i  ),ADD_VALUES);
                };
                if (fabs(Rhs(2*i+1)) >= 1.e-15){
                    int dof_i = 2 * connec(i) + 1;
                    ierr = VecSetValues(b,1,&dof_i,&Rhs(2*i+1),ADD_VALUES);
                };
                if (fabs(Rhs(12+i)) >= 1.e-15){
                    int dof_i = 2 * numNodesCoarse + connec(i);
                    ierr = VecSetValues(b,1,&dof_i,&Rhs(12+i),ADD_VALUES);
                };
            };
        };//Domain decomposition Coarse
    };//Elements Coarse
          
    //Fine mesh
    for (int jel = 0; jel < numElemFine; jel++){   
        
        if (domDecompFine.first[jel] == rank) {
                              
            //Compute Element matrix
            elementsFine_[jel] -> getTransientNavierStokesInitial();
          
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
                        int dof_i = 3*numNodesCoarse + 2 * connec(i);
                        int dof_j = 3*numNodesCoarse + 2 * connec(j);
                        ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&Ajac(2*i  ,2*j  ),ADD_VALUES);
                    };
                    if (fabs(Ajac(2*i+1,2*j  )) >= 1.e-15){
                        int dof_i = 3*numNodesCoarse + 2 * connec(i) +1;
                        int dof_j = 3*numNodesCoarse + 2 * connec(j);
                        ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&Ajac(2*i+1,2*j  ),ADD_VALUES);
                    };
                    if (fabs(Ajac(2*i  ,2*j+1)) >= 1.e-15){
                        int dof_i = 3*numNodesCoarse + 2 * connec(i);
                        int dof_j = 3*numNodesCoarse + 2 * connec(j) +1;
                        ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&Ajac(2*i  ,2*j+1),ADD_VALUES);
                    };
                    if (fabs(Ajac(2*i+1,2*j+1)) >= 1.e-15){
                        int dof_i = 3*numNodesCoarse + 2 * connec(i) +1;
                        int dof_j = 3*numNodesCoarse + 2 * connec(j) +1;
                        ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&Ajac(2*i+1,2*j+1),ADD_VALUES);
                    };
                    //Matrix Q and Qt
                    if (fabs(Ajac(2*i  ,12+j)) >= 1.e-15){
                        int dof_i = 3 * numNodesCoarse + 2 * connec(i);
                        int dof_j = 3 * numNodesCoarse + 2 * numNodesFine + connec(j);
                        ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&Ajac(2*i  ,12+j),ADD_VALUES);
                    };
                    if (fabs(Ajac(12+j,2*i  )) >= 1.e-15){
                        int dof_i = 3 * numNodesCoarse + 2 * connec(i);
                        int dof_j = 3 * numNodesCoarse + 2 * numNodesFine + connec(j);
                        ierr = MatSetValues(A,1,&dof_j,1,&dof_i,&Ajac(12+j,2*i  ),ADD_VALUES);
                    };
                    if (fabs(Ajac(2*i+1,12+j)) >= 1.e-15){
                        int dof_i = 3 * numNodesCoarse + 2 * connec(i) + 1;
                        int dof_j = 3 * numNodesCoarse + 2 * numNodesFine + connec(j);
                        ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&Ajac(2*i+1,12+j),ADD_VALUES);
                    };           
                    if (fabs(Ajac(12+j,2*i+1)) >= 1.e-15){
                        int dof_i = 3 * numNodesCoarse + 2 * connec(i) + 1;
                        int dof_j = 3 * numNodesCoarse + 2 * numNodesFine + connec(j);
                        ierr = MatSetValues(A,1,&dof_j,1,&dof_i,&Ajac(12+j,2*i+1),ADD_VALUES);
                    };
                    if (fabs(Ajac(12+i,12+j)) >= 1.e-15){
                        int dof_i = 3 * numNodesCoarse + 2 * numNodesFine + connec(i);
                        int dof_j = 3 * numNodesCoarse + 2 * numNodesFine + connec(j);
                        ierr = MatSetValues(A,1,&dof_j,1,&dof_i,&Ajac(12+i,12+j),ADD_VALUES);
                    };
                };     
                ///Rhs vector
                if (fabs(Rhs(2*i  )) >= 1.e-15){
                    int dof_i = 3 * numNodesCoarse + 2 * connec(i);
                    ierr = VecSetValues(b,1,&dof_i,&Rhs(2*i  ),ADD_VALUES);
                };
                if (fabs(Rhs(2*i+1)) >= 1.e-15){
                    int dof_i = 3 * numNodesCoarse + 2 * connec(i) + 1;
                    ierr = VecSetValues(b,1,&dof_i,&Rhs(2*i+1),ADD_VALUES);
                };
                if (fabs(Rhs(12+i)) >= 1.e-15){
                    int dof_i = 3 * numNodesCoarse + 2 * numNodesFine + connec(i);
                    ierr = VecSetValues(b,1,&dof_i,&Rhs(12+i),ADD_VALUES);
                };
            }; 
        };//Domain decomposition Fine             
    };//Elements Fine
                              
    //Lagrange Multipliers
    for (int l=0; l< numElemGlueZoneFine; l++){
                          
        int jel = elementsGlueZoneFine_[l];
                          
        if (domDecompFine.first[jel] == rank) {
          
            typename Elements::LocalMatrix Ajac, AjacAnt, AStab, ArlequinM, ArlequinM2;
            typename Elements::LocalVector Rhs, RhsStab, RhsArlequin,rhsLagMult;
            typename Elements::Connectivity connec, connecC, connecL;
              
            connec = elementsFine_[jel] -> getConnectivity();
            connecL = glueZoneFine_[l] -> getConnectivity();

            // FINE MESH
            
            //Computes element matrix
            elementsFine_[jel] -> getLagrangeMultipliersSameMesh();
              
            //Gets element matrice and rhs vectors
            // -L1
            Ajac = - elementsFine_[jel] -> getLagrMultMatrix();
            // +L1 * Lambda
            rhsLagMult = -elementsFine_[jel] -> getRhsLagrangeMultipliers(-Ajac);
            AStab = elementsFine_[jel] -> getJacNRMatrix();
            RhsStab = elementsFine_[jel] -> getRhsVector();
            ArlequinM = elementsFine_[jel] -> getArlequinStabilizationMatrix();
            RhsArlequin = elementsFine_[jel] -> getArlequinStabilizationVector();
              
            ArlequinM2 = elementsFine_[jel] -> getArlequinStabilizationMatrix2();
                
            // -L1t
            Ajac = trans(Ajac);
              
            // AStab.clear();
            // RhsStab.clear();
            // ArlequinM.clear();
            // RhsArlequin.clear();

            // elementsFine_[jel] -> getLMStabilizationSameMesh();

            // AStab = - elementsFine_[jel] -> getJacNRMatrix();
            // RhsStab = - elementsFine_[jel] -> getRhsVector();

            //AjacAnt = trans(AjacAnt);
              
              
            // U_.clear();
            // rhsLagMult.clear();
            // lagStab.clear();
            // for (int i = 0; i < 6; i++){
            //     U_(2*i  ) = nodesFine_[connec(i)] -> getLagrangeMultiplier(0);
            //     U_(2*i+1) = nodesFine_[connec(i)] -> getLagrangeMultiplier(1);
            // };
              
            // noalias(rhsLagMult) = -prod((Ajac),U_);
            // noalias(lagStab) = -prod(trans(AjacAnt),U_);
            // L1t * u1
            Rhs = -elementsFine_[jel] -> getRhsVelocities(-Ajac);

            std::pair<Elements::LocalVector, Elements::LocalMatrix>  lagMult;

            lagMult.first.clear();
            lagMult.second.clear();
                   
            lagMult = elementsFine_[jel] -> getBoundaryConditionsVelocity(Rhs,Ajac);

            Ajac = lagMult.second;
            Rhs = lagMult.first;

            lagMult.first.clear();
            lagMult.second.clear();
                   
            lagMult = elementsFine_[jel] -> getBoundaryConditionsLagMult(rhsLagMult,Ajac);

            Ajac = lagMult.second;
            rhsLagMult = lagMult.first;

            // for (int i=0; i<6; i++){
            //     if(lagMult.second(12+i,12+i) > 0){
            //         rhsLagMult(2*i  ) = 0.;
            //         rhsLagMult(2*i+1) = 0.;
            //         double one = 1.;
            //         int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
            //         ierr = MatSetValues(A,1,&d_i,1,&d_i,&one,INSERT_VALUES);
            //         d_i++;
            //         ierr = MatSetValues(A,1,&d_i,1,&d_i,&one,INSERT_VALUES);
            //     }
            // }
            double integ = alpha_f * gamma * dTime;
            //Disperse local contributions into the global matrix
            for (int i=0; i<6; i++){
                for (int j=0; j<6; j++){   
                    //COUPLING OPERATOR
                    if (fabs(Ajac(2*i  ,2*j  )) >= 1.e-15){
                        int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                        int d_j = 3*numNodesCoarse + 2*connec(j);
                        double value = Ajac(2*i  ,2*j  )*integ;
                        ierr = MatSetValues(A,1,&d_i,1,&d_j,&value,ADD_VALUES);
                        ierr = MatSetValues(A,1,&d_j,1,&d_i,&Ajac(2*i  ,2*j  ),ADD_VALUES);
                    };
                    if (fabs(Ajac(2*i+1,2*j  )) >= 1.e-15){
                        int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i) + 1;
                        int d_j = 3*numNodesCoarse + 2*connec(j);
                        double value = Ajac(2*i+1,2*j  )*integ;
                        ierr = MatSetValues(A,1,&d_i,1,&d_j,&value,ADD_VALUES);
                        ierr = MatSetValues(A,1,&d_j,1,&d_i,&Ajac(2*i+1,2*j  ),ADD_VALUES);
                    };
                    if (fabs(Ajac(2*i+1,2*j+1)) >= 1.e-15){
                        int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i) + 1;
                        int d_j = 3*numNodesCoarse + 2*connec(j) + 1;
                        double value = Ajac(2*i+1,2*j+1)*integ;
                        ierr = MatSetValues(A,1,&d_i,1,&d_j,&value,ADD_VALUES);
                        ierr = MatSetValues(A,1,&d_j,1,&d_i,&Ajac(2*i+1,2*j+1),ADD_VALUES);
                    };
                    if (fabs(Ajac(2*i  ,2*j+1)) >= 1.e-15){
                        int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                        int d_j = 3 * numNodesCoarse + 2*connec(j) + 1;
                        double value = Ajac(2*i  ,2*j+1)*integ;
                        ierr = MatSetValues(A,1,&d_i,1,&d_j,&value,ADD_VALUES);
                        ierr = MatSetValues(A,1,&d_j,1,&d_i,&Ajac(2*i  ,2*j+1),ADD_VALUES);
                    };
                    if (fabs(Ajac(12+i,12+j)) >= 1.e-15){
                        int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                        int d_j = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(j);                      
                        //std::cout << "AQUI " << Ajac(12+i,12+j)<< " " << i << " " << j << " " << d_i << " " << d_j << std::endl;
                        // ierr = MatSetValues(A,1,&d_i,1,&d_j,&Ajac(12+i,12+j),ADD_VALUES);
                        // d_i++; d_j++;
                        // ierr = MatSetValues(A,1,&d_i,1,&d_j,&Ajac(12+i,12+j),ADD_VALUES);
                    };
                    //SUPG STABILIZATION
                    if (fabs(AStab(2*i  ,2*j  )) >= 1.e-15){
                        int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                        int d_j = 3*numNodesCoarse + 2*connec(j);
                        double value = AStab(2*i  ,2*j  )*integ;
                        ierr = MatSetValues(A,1,&d_j,1,&d_i,&value,ADD_VALUES);
                    };
                    if (fabs(AStab(2*i+1,2*j+1)) >= 1.e-15){
                        int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i) + 1;
                        int d_j = 3*numNodesCoarse + 2*connec(j) + 1;
                        double value = AStab(2*i+1,2*j+1)*integ;
                        ierr = MatSetValues(A,1,&d_j,1,&d_i,&value,ADD_VALUES);
                    };
                    //PSPG STABILIZATION
                    if (fabs(AStab(2*i  ,12+j)) >= 1.e-15){
                        int dof_i = 3*numNodesCoarse + 2*numNodesFine + connec(i);
                        int dof_j = 3*numNodesCoarse + 3*numNodesFine + connecL(j);
                        ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&AStab(2*i  ,12+j),ADD_VALUES);
                    };
                    if (fabs(AStab(2*i+1,12+j)) >= 1.e-15){
                        int dof_i = 3*numNodesCoarse + 2*numNodesFine + connec(i);
                        int dof_j = 3*numNodesCoarse + 3*numNodesFine + connecL(j)+1;
                        ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&AStab(2*i+1,12+j),ADD_VALUES);
                    };

                    //ARLEQUIN STABILIZATION
                    // if (fabs(ArlequinM(2*i  ,2*j+1)) >= 1.e-15){
                    //     int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                    //     int d_j = 3*numNodesCoarse + 2*connec(j) + 1;
                    //     ierr = MatSetValues(A,1,&d_i,1,&d_j,
                    //                         &ArlequinM(2*i  ,2*j+1),ADD_VALUES);
                    // };
                    if (fabs(ArlequinM2(2*i  ,2*j  )) >= 1.e-15){
                        int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                        int d_j = 3*numNodesCoarse + 2*connec(j);
                        ierr = MatSetValues(A,1,&d_i,1,&d_j,&ArlequinM2(2*i  ,2*j  ),ADD_VALUES);
                    };
                    // if (fabs(ArlequinM(2*i+1,2*j  )) >= 1.e-15){
                    //     int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i) + 1;
                    //     int d_j = 3*numNodesCoarse + 2*connec(j);
                    //     ierr = MatSetValues(A,1,&d_i,1,&d_j,
                    //                         &ArlequinM(2*i+1,2*j  ),ADD_VALUES);
                    // };
                    if (fabs(ArlequinM2(2*i+1,2*j+1)) >= 1.e-15){
                        int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i) + 1;
                        int d_j = 3*numNodesCoarse + 2*connec(j) + 1;
                        ierr = MatSetValues(A,1,&d_i,1,&d_j,&ArlequinM2(2*i+1,2*j+1),ADD_VALUES);
                    };
                    if (fabs(ArlequinM2(12+i,2*j  )) >= 1.e-15){
                        int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                        int dof_j = 3*numNodesCoarse + 2*numNodesFine + connec(j);
                        ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinM2(12+i,2*j  ),ADD_VALUES);
                    };
                    if (fabs(ArlequinM2(12+i,2*j+1)) >= 1.e-15){
                        int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i) + 1;
                        int dof_j = 3*numNodesCoarse + 2*numNodesFine + connec(j);
                        ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinM2(12+i,2*j+1),ADD_VALUES);
                    };
                    if (fabs(ArlequinM(2*i  ,2*j  )) >= 1.e-15){
                        int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                        int dof_j = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(j);
                        ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinM(2*i  ,2*j  ),ADD_VALUES);
                        dof_i++; dof_j++;
                        ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinM(2*i+1,2*j+1),ADD_VALUES);
                    };                            
                };
                //RHS VECTOR
                //COUPLING OPERATOR
                int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                ierr = VecSetValues(b,1,&dof_i,&Rhs(2*i  ),ADD_VALUES);

                dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i) + 1;
                ierr = VecSetValues(b,1,&dof_i,&Rhs(2*i+1),ADD_VALUES);

                dof_i = 3*numNodesCoarse + 2*connec(i);
                ierr = VecSetValues(b,1,&dof_i,&rhsLagMult(2*i  ),ADD_VALUES);

                dof_i = 3*numNodesCoarse + 2*connec(i) + 1;
                ierr = VecSetValues(b,1,&dof_i,&rhsLagMult(2*i+1),ADD_VALUES);

                //SUPG STABILIZATION
                dof_i = 3*numNodesCoarse + 2*connec(i);                       
                ierr = VecSetValues(b,1,&dof_i,&RhsStab(2*i  ),ADD_VALUES);

                dof_i = 3*numNodesCoarse + 2*connec(i) + 1;
                ierr = VecSetValues(b,1,&dof_i,&RhsStab(2*i+1),ADD_VALUES);

                //PSPG STABILIZATION
                dof_i = 3*numNodesCoarse + 2*numNodesFine + connec(i);
                ierr = VecSetValues(b,1,&dof_i,&RhsStab(12+i),ADD_VALUES);

                // ///ARLEQUIN STABILIZATION
                dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                ierr = VecSetValues(b,1,&dof_i,&RhsArlequin(2*i  ),ADD_VALUES);
                dof_i++;
                ierr = VecSetValues(b,1,&dof_i,&RhsArlequin(2*i+1),ADD_VALUES);
            };      
              
            //COAESE MESH
                                  
            //Counts number of coarse mesh intersecting the 
            //fine element
            int numberIntPoints = elementsFine_[jel] -> getNumberOfIntegrationPoints();
            int aux;
              
            std::vector<int> ele, diffElem;
            ele.clear();
            diffElem.clear();

            ele.reserve(3);
            for (int i=0; i<numberIntPoints; i++){
                aux = elementsFine_[jel] -> getIntegPointCorrespondenceElement(i);
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
                double pspg = 0.;//elementsCoarse_[iElemCoarse] -> getPSPG();
                ublas::bounded_vector<double, 6> press_, velX_, velY_;

                connecC = elementsCoarse_[iElemCoarse] -> getConnectivity();

                for (int k = 0; k < 6; k++){
                    press_(k) = nodesCoarse_[connecC(k)] -> getPressure();
                    velX_(k) = nodesCoarse_[connecC(k)] -> getVelocity(0);
                    velY_(k) = nodesCoarse_[connecC(k)] -> getVelocity(1);
                }
                  
                elementsFine_[jel] -> getLagrangeMultipliersDifferentMesh(iElemCoarse,pspg,press_,velX_,velY_);
                  
                //Computes element matrix
                // L0
                Ajac = elementsFine_[jel] -> getLagrMultMatrix();
                // L0t
                Ajac = trans(Ajac);
                AStab = elementsFine_[jel] -> getJacNRMatrix();
                AStab = trans(AStab);
                RhsStab = elementsFine_[jel] -> getRhsVector();
                ArlequinM = elementsFine_[jel] -> getArlequinStabilizationMatrix();
                ArlequinM = trans(ArlequinM);
                RhsArlequin = elementsFine_[jel] -> getArlequinStabilizationVector();
                  
                ArlequinM2 = elementsFine_[jel] -> getArlequinStabilizationMatrix2();
                ArlequinM2 = trans(ArlequinM2);

                // -L0 * lambda
                rhsLagMult =  elementsFine_[jel] -> getRhsLagrangeMultipliers(trans(Ajac));


                // AStab.clear();
                // // RhsStab.clear();
                // ArlequinM.clear();
                // RhsArlequin.clear();

                // U_.clear();
                // for (int i = 0; i < 6; i++){
                //     U_(2*i  ) = nodesFine_[connec(i)] -> getLagrangeMultiplier(0);
                //     U_(2*i+1) = nodesFine_[connec(i)] -> getLagrangeMultiplier(1);
                // };
                  
                // noalias(rhsLagMult) = -prod(trans(Ajac),U_);
                // noalias(lagStab) = -prod(trans(AjacAnt),U_);
                  
                Rhs = elementsCoarse_[iElemCoarse] -> getRhsVelocities(Ajac);

                std::pair<Elements::LocalVector, Elements::LocalMatrix>  lagMult;
              
                // lagMult.first.clear();
                // lagMult.second.clear();
            
                // lagMult = elementsCoarse_[iElemCoarse] -> getBoundaryConditionsVelocity(rhsLagMult,(Ajac));

                // Ajac = lagMult.second;
                // rhsLagMult = lagMult.first;

                // lagMult.first.clear();
                // lagMult.second.clear();
                   
                // lagMult = elementsFine_[jel] -> getBoundaryConditionsLagMult(Rhs,(Ajac));

                // Ajac = (lagMult.second);
                // Rhs = lagMult.first;

                // if (jel == 94){
                //     std::cout << "MATRIZ " << std::endl;
                //     for (int i=0; i<18; i++){
                //         for (int j=0; j<18; j++){
                //             std::cout << Ajac(i,j) << " ";
                //         }
                //         std::cout << std::endl;
                //     }
                //     std::cout << "VETOR " << std::endl;
                //     for (int i=0; i<18; i++){
                //         std::cout << Rhs(i) << " " << rhsLagMult(i) << std::endl;
                //     }
                // }

                rhsLagMult *= 1000.e0;
                //Ajac *= 1000.;

                for (int i = 0; i < 12; i++){
                    RhsStab(i) *= 1000.e0;
                    for (int j = 0; j < 12; j++){
                        AStab(i,j)*=1000.e0;
                        Ajac(i,j)*=1000.e0;
                        ArlequinM2(i,j)*=1000.e0;
                    }
                }
                // std::pair<Elements::LocalVector, 
                //           Elements::LocalMatrix>  lagMult;
                  
                // lagMult.first.clear();
                // lagMult.second.clear();
                   
                // lagMult = elementsCoarse_[iElemCoarse] -> 
                //     getRhsVectorAndBoundaryConditions(Ajac);

                // noalias(rhsLagMult) = -prod((lagMult.second),U_);
                // for (int i=0; i<6; i++){
                //     if(lagMult.second(12+i,12+i) > 0){
                //         rhsLagMult(2*i  ) = 0.;
                //         rhsLagMult(2*i+1) = 0.;
                //         double one = 1.;
                //         int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                //         ierr = MatSetValues(A,1,&d_i,1,&d_i,&one,INSERT_VALUES);
                //         d_i++;
                //         ierr = MatSetValues(A,1,&d_i,1,&d_i,&one,INSERT_VALUES);
                //     }
                // }

                double integ = alpha_f * gamma * dTime;

                //Disperse local contribution into the global matrix
                for (int i=0; i<6; i++){
                    for (int j=0; j<6; j++){
                        //COUPLING OPERATOR
                        if (fabs(Ajac(2*i  ,2*j  )) >= 1.e-15){
                            int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                            int dof_j = 2*connecC(j);
                            double value = Ajac(2*i  ,2*j  )*integ;
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&value,ADD_VALUES);
                            ierr = MatSetValues(A,1,&dof_j,1,&dof_i,&Ajac(2*i  ,2*j  ),ADD_VALUES);
                        }
                        if (fabs(Ajac(2*i+1,2*j  )) >= 1.e-15){
                            int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i) + 1;
                            int dof_j = 2*connecC(j);
                            double value = Ajac(2*i+1,2*j  )*integ;
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&value,ADD_VALUES);
                            ierr = MatSetValues(A,1,&dof_j,1,&dof_i,&Ajac(2*i+1,2*j  ),ADD_VALUES);
                        };
                        if (fabs(Ajac(2*i+1,2*j+1)) >= 1.e-15){
                            int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i) + 1;
                            int dof_j = 2*connecC(j) + 1;
                            double value = Ajac(2*i+1,2*j+1)*integ;
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&value,ADD_VALUES);
                            ierr = MatSetValues(A,1,&dof_j,1,&dof_i,&Ajac(2*i+1,2*j+1),ADD_VALUES);
                        };
                        if (fabs(Ajac(2*i  ,2*j+1)) >= 1.e-15){
                            int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                            int dof_j = 2*connecC(j) + 1;
                            double value = Ajac(2*i  ,2*j+1)*integ;
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&value,ADD_VALUES);
                            ierr = MatSetValues(A,1,&dof_j,1,&dof_i,&Ajac(2*i  ,2*j+1),ADD_VALUES);
                        };
                        if (fabs(Ajac(12+i,12+j)) >= 1.e-15){
                            int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                            int dof_j = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(j);
                            // ierr = MatSetValues(A,1,&dof_i,1,&dof_j,
                            //                 &lagMult.second(12+i,12+j),ADD_VALUES);
                            dof_i++; dof_j++;
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&Ajac(12+i,12+j),ADD_VALUES);
                        };
                        //SUPG STABILIZATION
                        if (fabs(AStab(2*i  ,2*j  )) >= 1.e-15){
                            int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                            int dof_j = 2*connecC(j);
                            double value = AStab(2*i  ,2*j  )*integ;
                            ierr = MatSetValues(A,1,&dof_j,1,&dof_i,&value,ADD_VALUES);
                        };
                        if (fabs(AStab(2*i+1,2*j+1)) >= 1.e-15){
                            int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i) + 1;
                            int dof_j = 2*connecC(j) + 1;
                            double value = AStab(2*i+1,2*j+1) * integ;
                            ierr = MatSetValues(A,1,&dof_j,1,&dof_i,&value,ADD_VALUES);
                        };
                        //PSPG STABILIZATION
                        if (fabs(AStab(2*i  ,12+j)) >= 1.e-15){
                            int dof_i = 2*numNodesCoarse + connecC(i);
                            int dof_j = 3*numNodesCoarse + 2*numNodesFine + connecL(j);
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&AStab(2*i  ,12+j),ADD_VALUES);
                        };
                        if (fabs(AStab(2*i+1,12+j)) >= 1.e-15){
                            int dof_i = 2*numNodesCoarse + connecC(i);
                            int dof_j = 3*numNodesCoarse + 2*numNodesFine + connecL(j)+1;
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&AStab(2*i+1,12+j),ADD_VALUES);
                        };
                        ///ARLEQUIN STABILIZATION
                        // if (fabs(ArlequinM(2*i  ,2*j+1)) >= 1.e-15){
                        //     int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                        //     int d_j = 2*connecC(j) + 1;
                        //     ierr = MatSetValues(A,1,&d_i,1,&d_j,
                        //                         &ArlequinM(2*i  ,2*j+1),ADD_VALUES);
                        // };
                        if (fabs(ArlequinM2(2*i  ,2*j  )) >= 1.e-15){
                            int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                            int d_j = 2*connecC(j);
                            ierr = MatSetValues(A,1,&d_i,1,&d_j,&ArlequinM2(2*i  ,2*j  ),ADD_VALUES);
                        };
                        // if (fabs(ArlequinM(2*i+1,2*j  )) >= 1.e-15){
                        //     int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i) + 1;
                        //     int d_j = 2*connecC(j);
                        //     ierr = MatSetValues(A,1,&d_i,1,&d_j,
                        //                         &ArlequinM(2*i+1,2*j  ),ADD_VALUES);
                        // };
                        if (fabs(ArlequinM2(2*i+1,2*j+1)) >= 1.e-15){
                            int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i) + 1;
                            int d_j = 2*connecC(j) + 1;
                            ierr = MatSetValues(A,1,&d_i,1,&d_j,&ArlequinM2(2*i+1,2*j+1),ADD_VALUES);
                        };
                        if (fabs(ArlequinM2(12+i,2*j  )) >= 1.e-15){
                            int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                            int dof_j = 2*numNodesCoarse + connecC(j);
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinM2(12+i,2*j  ),ADD_VALUES);
                        };
                        if (fabs(ArlequinM2(12+i,2*j+1)) >= 1.e-15){
                            int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i) + 1;
                            int dof_j = 2*numNodesCoarse + connecC(j);
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinM2(12+i,2*j+1),ADD_VALUES);
                        };
                        if (fabs(ArlequinM(2*i  ,2*j  )) >= 1.e-15){
                            int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                            int dof_j = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(j);
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinM(2*i  ,2*j  ),ADD_VALUES);
                            dof_i++; dof_j++;
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinM(2*i+1,2*j+1),ADD_VALUES);
                        };
                    };
                    //RHS VECTOR
                    //COUPLING OPERATOR
                    int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                    ierr = VecSetValues(b,1,&d_i,&Rhs(2*i  ),ADD_VALUES);

                    d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i) + 1;
                    ierr = VecSetValues(b,1,&d_i,&Rhs(2*i+1),ADD_VALUES);

                    int dof_i = 2*connecC(i);
                    ierr = VecSetValues(b,1,&dof_i,&rhsLagMult(2*i  ),ADD_VALUES);

                    dof_i = 2*connecC(i) + 1;
                    ierr = VecSetValues(b,1,&dof_i,&rhsLagMult(2*i+1),ADD_VALUES);

                    //SUPG STABILIZATION
                    dof_i = 2*connecC(i);
                    ierr = VecSetValues(b,1,&dof_i,&RhsStab(2*i  ),ADD_VALUES);

                    dof_i = 2*connecC(i) + 1;
                    ierr = VecSetValues(b,1,&dof_i,&RhsStab(2*i+1),ADD_VALUES);

                    //PSPG STABILIZATION
                    dof_i = 2*numNodesCoarse + connecC(i);
                    ierr = VecSetValues(b,1,&dof_i,&RhsStab(12+i),ADD_VALUES);

                    //ARLEQUIN STABILIZATION
                    dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                    ierr = VecSetValues(b,1,&dof_i,&RhsArlequin(2*i  ),ADD_VALUES);
                    dof_i++;
                    ierr = VecSetValues(b,1,&dof_i,&RhsArlequin(2*i+1),ADD_VALUES);
                };                                 
            }; //Number of intersections
        }; // if element belongs to the glue zone
    }; // Glue zone
          
                    
    //Assemble matrices and vectors
    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
                      
    ierr = VecAssemblyBegin(b);
    ierr = VecAssemblyEnd(b);
                                
    //Create KSP context to solve the linear system
    ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);
                      
    ierr = KSPSetOperators(ksp,A,A);
                      
                      
#if defined(PETSC_HAVE_MUMPS)
    ierr = KSPSetType(ksp,KSPPREONLY);
    ierr = KSPGetPC(ksp,&pc);
    ierr = PCSetType(pc, PCLU);
                      
    ierr = PCFactorSetMatSolverType(pc,MATSOLVERMUMPS);
    PCFactorSetUpMatSolverType(pc);
    PCFactorGetMatrix(pc,&F);
                      
    PetscInt ival,icntl;
    icntl = 14; ival = 80;
    MatMumpsSetIcntl(F,icntl,ival);
    // icntl = 11; ival = 1;
    // MatMumpsSetIcntl(F,11,1);
          
    //MatMumpsSetIcntl(F,21,0);       
#endif
    ierr = KSPSetFromOptions(ksp);
    ierr = KSPSetUp(ksp);
                      
#if defined(PETSC_HAVE_MUMPS)
    PetscInt  info1,info2,icntl14;
                      
    MatMumpsGetInfo(F,1,&info1);
    MatMumpsGetInfo(F,2,&info2);
          
    // PetscReal info5,info6,info7,info8,info9,info10,info11;
    // MatMumpsGetRinfog(F,5,&info5);
    // MatMumpsGetRinfog(F,6,&info6);
    // MatMumpsGetRinfog(F,7,&info7);
    // MatMumpsGetRinfog(F,8,&info8);
    // MatMumpsGetRinfog(F,9,&info9);
    // MatMumpsGetRinfog(F,10,&info10);
    // MatMumpsGetRinfog(F,11,&info11);

    // PetscInt info21,info32;
    // MatMumpsGetIcntl(F,32,&info32);
    // MatMumpsGetIcntl(F,21,&info21);
          
                      
    // if(rank==0) std::cout << "ICNTL = " << info21 << " " << info32 << std::endl;

    // if(rank==0) std::cout << "INFOG = " << info5 << " " << info6 << " " << info7 << " " << info8 << " " << info9 << " " << info10 << " " << info11 << std::endl;
    MatMumpsGetIcntl(F,14,&icntl14);    
    if((rank == 0) && (info1 != 0)) std::cout << " INFO(1) = " << info1 << " " << info2 << " " 
                                                << icntl14 << std::endl;
#endif
                      
    // ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
      
    ierr = KSPSolve(ksp,b,u);
      
    ierr = KSPGetTotalIterations(ksp, &iterations);
      
    //if (rank == 0)std::cout << "GMRES Iterations = " << iterations << std::endl;
  
    //ierr = VecView(u,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
      
    //Gathers the solution vector to the master process
    ierr = VecScatterCreateToAll(u, &ctx, &All);  
    ierr = VecScatterBegin(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);
    ierr = VecScatterEnd(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);
    ierr = VecScatterDestroy(&ctx);
      
    //Updates nodal values
    double u_[2];
    double normU = 0.;
    double normP = 0.;
    double normL = 0.;
    double normT = 0.;
    double p_;
    Ione = 1;

    for (int i = 0; i < numNodesCoarse; ++i){
        Ii = 2 * i;
        ierr = VecGetValues(All, Ione, &Ii, &val);
        if (nodesCoarse_[i] -> getDistFunction() > -1.2) val *= 1000.e0;
        u_[0] = nodesCoarse_[i] -> getVelocity(0);
        normU += val*val;
        nodesCoarse_[i] -> setAccelerationComponent(0,val);
        nodesCoarse_[i] -> setPreviousVelocityComponent(0,u_[0]-val*gamma*dTime);
          
        Ii = 2 * i + 1;
        ierr = VecGetValues(All, Ione, &Ii, &val);
        if (nodesCoarse_[i] -> getDistFunction() > -1.2) val *= 1000.e0;
        u_[1] = nodesCoarse_[i] -> getVelocity(1);
        normU += val*val;
        nodesCoarse_[i] -> setAccelerationComponent(1,val);
        nodesCoarse_[i] -> setPreviousVelocityComponent(1,u_[1]-val*gamma*dTime);
    };     
    for (int i = 0; i<numNodesCoarse; i++){
        Ii = 2 * numNodesCoarse + i;
        ierr = VecGetValues(All,Ione,&Ii,&val);
        p_ = val;
        normP += val*val;
        nodesCoarse_[i] -> incrementPressure(p_);
    };
      
    for (int i = 0; i < numNodesFine; ++i){
        Ii = 3 * numNodesCoarse + 2 * i;
        ierr = VecGetValues(All, Ione, &Ii, &val);
        u_[0] = nodesFine_[i] -> getVelocity(0);
        normU += val*val;
        nodesFine_[i] -> setAccelerationComponent(0,val);
        nodesFine_[i] -> setPreviousVelocityComponent(0,u_[0]-val*gamma*dTime);
          
        Ii = 3 * numNodesCoarse + 2 * i + 1;
        ierr = VecGetValues(All, Ione, &Ii, &val);
        u_[1] = nodesFine_[i] -> getVelocity(1);
        normU += val*val;
        nodesFine_[i] -> setAccelerationComponent(1,val);
        nodesFine_[i] -> setPreviousVelocityComponent(1,u_[1]-val*gamma*dTime);
    };
      
    for (int i = 0; i<numNodesFine; i++){
        Ii = 3 * numNodesCoarse + 2 * numNodesFine + i;
        ierr = VecGetValues(All,Ione,&Ii,&val);
        p_ = val;
        normP += val*val;
        nodesFine_[i] -> incrementPressure(p_);
    };
      
    for (int i = 0; i < numNodesGlueZoneFine; ++i){
        Ii = 3 * numNodesCoarse + 3 * numNodesFine + 2 * i;
        ierr = VecGetValues(All, Ione, &Ii, &val);
        u_[0] = val;
        nodesFine_[nodesGlueZoneFine_[i]] -> incrementLagrangeMultiplier(0,u_[0]);
        normL += val*val;

        Ii = 3 * numNodesCoarse + 3 * numNodesFine + 2 * i + 1;
        ierr = VecGetValues(All, Ione, &Ii, &val);
        u_[1] = val;
        nodesFine_[nodesGlueZoneFine_[i]] -> incrementLagrangeMultiplier(1,u_[1]);
        normL += val*val;
        // std::cout << "LAG M " << u_[0] << " " << u_[1] << std::endl;
    };
                                        
    boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
      
    ierr = KSPDestroy(&ksp); 
    ierr = VecDestroy(&b); 
    ierr = VecDestroy(&u); 
    ierr = VecDestroy(&All); 
    ierr = MatDestroy(&A); 
    // ierr = MatDestroy(&F); CHKERRQ(ierr);
    //ierr = MatDestroy(&C); CHKERRQ(ierr);
    
    //Printing results
    if (rank == 0) printResults(100);
        
    return;

};






//------------------------------------------------------------------------------
//----------------COMPUTE ARLEQUIN COUPLED NAVIER-STOKES PROBLEM----------------
//------------------------------------------------------------------------------
template<>
int Arlequin<2>::solveArlequinProblem(int iterNumber, double tolerance,
                                      int problem_type, int time_dependency){

    Mat               A, F;
    Vec               b, u, All;
    PetscErrorCode    ierr;
    PetscInt          Istart, Iend, Ii, Ione, iterations;
    KSP               ksp; 
    PC                pc;
    VecScatter        ctx;
    PetscScalar       val;
    PetscLogDouble bytes = 0;


    std::ofstream dragLift;
    dragLift.open("dragLift.dat", std::ofstream::out | std::ofstream::app);
    if (rank == 0) {
        dragLift << "Time   Pressure Drag   Pressure Lift " 
                 << "Friction Drag  Friction Lift Drag    Lift " 
                 << std::endl;
    };

    if ((problem_type < 1) || (problem_type > 2)){
        std::cout << "WRONG PROBLEM TYPE." << std::endl;
        return 0;
    };

    if (time_dependency == 0) numTimeSteps = 1;
    
    iTimeStep = 0.;
   
    // //Computes the Weight function for all the finite elements
    setWeightFunction(16.);

    //Computes the Nodal correspondence between fine nodes and coarse elements
    setNodalCorrespondenceFine();
    // if (rank == 0) printResults(100);
    // Computes the system size
    int sysSize = 3 * numNodesCoarse +  
        + 3 * numNodesFine + 2 * numNodesGlueZoneFine;

    double integScheme = fineModel.integScheme;

    alpha_f = 1. / (1. + integScheme);
    alpha_m = 0.5 * (3. - integScheme) / (1. + integScheme);
    gamma = 0.5 + alpha_m - alpha_f;

    // initialAcceleration();

    for (iTimeStep = 0; iTimeStep < numTimeSteps; iTimeStep++){


        // if (iTimeStep <= 10){
        //     iterNumber = 3;
        // } else {
        //     iterNumber = 3;
        // }

        // if (iTimeStep == 25){
        //     for (int i = 0; i < numElemFine; ++i){
        //         double inSch = 0.5;
        //         elementsFine_[i] -> setTimeIntegrationScheme(inSch);
        //     }
        //     for (int i=0; i<numElemCoarse; i++){
        //         double inSch = 0.5;
        //         elementsCoarse_[i] -> setTimeIntegrationScheme(inSch);
        //     }
        // }
           
        parametersCoarse -> setTimeInstant(iTimeStep);
        parametersFine -> setTimeInstant(iTimeStep);

        if (rank == 0) {std::cout << "----------------------------" 
                                  << " TIME STEP = "
                                  << iTimeStep 
                                  << " ---------------------------"
                                  << std::endl;}
        PetscMemoryGetCurrentUsage(&bytes);
        PetscPrintf(PETSC_COMM_WORLD,"Memory used %g M\n",bytes/(1024*1024));

        //Updates velocity and acceleration
        for (int i = 0; i < numNodesCoarse; i++){
            double accel[2], u[2], uprev[2];
            
            //Compute acceleration
            u[0] = nodesCoarse_[i] -> getVelocity(0);
            u[1] = nodesCoarse_[i] -> getVelocity(1);

            nodesCoarse_[i] -> setPreviousVelocity(u);
            
            accel[0] = nodesCoarse_[i] -> getAcceleration(0);
            accel[1] = nodesCoarse_[i] -> getAcceleration(1);
            
            nodesCoarse_[i] -> setPreviousAcceleration(accel);

            accel[0] *= (gamma - 1.) / gamma;
            accel[1] *= (gamma - 1.) / gamma;
            
            nodesCoarse_[i] -> setAcceleration(accel);            

        };

        for (int i = 0; i < numNodesFine; i++){
            double accel[2], u[2], uprev[2], lag[2];
            
            //Compute acceleration
            u[0] = nodesFine_[i] -> getVelocity(0);
            u[1] = nodesFine_[i] -> getVelocity(1);

            nodesFine_[i] -> setPreviousVelocity(u);
            
            accel[0] = nodesFine_[i] -> getAcceleration(0);
            accel[1] = nodesFine_[i] -> getAcceleration(1);
            
            nodesFine_[i] -> setPreviousAcceleration(accel);

            accel[0] *= (gamma - 1.) / gamma;
            accel[1] *= (gamma - 1.) / gamma;
            
            nodesFine_[i] -> setAcceleration(accel);
        };

        //STARTS NEWTON-RAPHSON
        for (int inewton = 0; inewton < iterNumber; inewton++){
            
            boost::posix_time::ptime t1 =
                boost::posix_time::microsec_clock::local_time();
            
            // Preallocates the matrix
            ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                                sysSize, sysSize, 300, NULL, 900, NULL, &A); 
            
            CHKERRQ(ierr);
            
            // Divides the matrix between the processes
            ierr = MatGetOwnershipRange(A, &Istart, &Iend);CHKERRQ(ierr);
            
            //Create PETSc vectors
            ierr = VecCreate(PETSC_COMM_WORLD, &b); CHKERRQ(ierr);
            ierr = VecSetSizes(b, PETSC_DECIDE, sysSize); CHKERRQ(ierr);
            ierr = VecSetFromOptions(b); CHKERRQ(ierr); 
            ierr = VecDuplicate(b, &u); CHKERRQ(ierr);
            ierr = VecDuplicate(b, &All); CHKERRQ(ierr);
                        
            for (int i=0; i<sysSize; i++){
                double val = 1.e-20;
                ierr = MatSetValues(A,1,&i,1,&i,&val,ADD_VALUES);
                
            }
                                

            //------------------------------------------------------------------
            //------------------BEGIN OF LINEAR SYSTEM ASSEMBLY-----------------
            //------------------------------------------------------------------

            //Coarse mesh
            
            for (int jel = 0; jel < numElemCoarse; jel++){   
                
                if (domDecompCoarse.first[jel] == rank) {
                    
                    //Compute Element matrix
                    std::pair<Elements::LocalMatrix, Elements::LocalVector> localMV;
                    localMV = elementsCoarse_[jel] -> getTransientNavierStokes();
        
                    typename Elements::Connectivity connec;
                    connec = elementsCoarse_[jel] -> getConnectivity();
                    
                    //Disperse local contributions into the global matrix
                    for (int i=0; i<6; i++){
                        for (int j=0; j<6; j++){
                            
                            if (fabs(localMV.first(2*i  ,2*j  )) >= 1.e-15){
                                int dof_i = 2 * connec(i);
                                int dof_j = 2 * connec(j);
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&localMV.first(2*i  ,2*j  ),ADD_VALUES);
                            };
                            if (fabs(localMV.first(2*i+1,2*j  )) >= 1.e-15){
                                int dof_i = 2 * connec(i) + 1;
                                int dof_j = 2 * connec(j);
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&localMV.first(2*i+1,2*j  ),ADD_VALUES);
                            };

                            if (fabs(localMV.first(2*i  ,2*j+1)) >= 1.e-15){
                                int dof_i = 2 * connec(i);
                                int dof_j = 2 * connec(j) + 1;
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&localMV.first(2*i  ,2*j+1),ADD_VALUES);
                            };

                            if (fabs(localMV.first(2*i+1,2*j+1)) >= 1.e-15){
                                int dof_i = 2 * connec(i) + 1;
                                int dof_j = 2 * connec(j) + 1;
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&localMV.first(2*i+1,2*j+1),ADD_VALUES);
                            };

                            //Matrix Q and Qt
                            if (fabs(localMV.first(2*i  ,12+j)) >= 1.e-15){
                                int dof_i = 2 * connec(i);
                                int dof_j = 2 * numNodesCoarse + connec(j);
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&localMV.first(2*i  ,12+j),ADD_VALUES);
                            };
                            
                            if (fabs(localMV.first(12+j,2*i  )) >= 1.e-15){
                                int dof_i = 2 * connec(i);
                                int dof_j = 2 * numNodesCoarse + connec(j);
                                ierr = MatSetValues(A,1,&dof_j,1,&dof_i,&localMV.first(12+j,2*i  ),ADD_VALUES);
                            };

                            if (fabs(localMV.first(2*i+1,12+j)) >= 1.e-15){
                                int dof_i = 2 * connec(i) + 1;
                                int dof_j = 2 * numNodesCoarse + connec(j);
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&localMV.first(2*i+1,12+j),ADD_VALUES);
                            };

                            if (fabs(localMV.first(12+j,2*i+1)) >= 1.e-15){
                                int dof_i = 2 * connec(i) + 1;
                                int dof_j = 2 * numNodesCoarse + connec(j);
                                ierr = MatSetValues(A,1,&dof_j,1,&dof_i,&localMV.first(12+j,2*i+1),ADD_VALUES);
                            };

                            if (fabs(localMV.first(12+i,12+j)) >= 1.e-15){
                                int dof_i = 2 * numNodesCoarse + connec(i);
                                int dof_j = 2 * numNodesCoarse + connec(j);
                                ierr = MatSetValues(A,1,&dof_j,1,&dof_i,&localMV.first(12+i,12+j),ADD_VALUES);
                            };
                        };                  
                        //Rhs vector
                        if (fabs(localMV.second(2*i  )) >= 1.e-15){
                            int dof_i = 2 * connec(i);
                            ierr = VecSetValues(b,1,&dof_i,&localMV.second(2*i  ),ADD_VALUES);
                        };

                        if (fabs(localMV.second(2*i+1)) >= 1.e-15){
                            int dof_i = 2 * connec(i) + 1;
                            ierr = VecSetValues(b,1,&dof_i,&localMV.second(2*i+1),ADD_VALUES);
                        };

                        if (fabs(localMV.second(12+i)) >= 1.e-15){
                            int dof_i = 2 * numNodesCoarse + connec(i);
                            ierr = VecSetValues(b,1,&dof_i,&localMV.second(12+i),ADD_VALUES);
                        };
                    };
                };
            };

            //Fine mesh
            for (int jel = 0; jel < numElemFine; jel++){   
                
                if (domDecompFine.first[jel] == rank) {
                    
                    //Compute Element matrix
                    std::pair<Elements::LocalMatrix, Elements::LocalVector> localMV;
                    localMV = elementsFine_[jel] -> getTransientNavierStokes();                    
                    
                    typename Elements::Connectivity connec;
                    connec = elementsFine_[jel] -> getConnectivity();
                    
                    //Disperse local contributions into the global matrix
                    for (int i=0; i<6; i++){
                        for (int j=0; j<6; j++){
                            if (fabs(localMV.first(2*i  ,2*j  )) >= 1.e-15){
                                int dof_i = 3*numNodesCoarse + 2 * connec(i);
                                int dof_j = 3*numNodesCoarse + 2 * connec(j);
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&localMV.first(2*i  ,2*j  ),ADD_VALUES);
                            };
                            
                            if (fabs(localMV.first(2*i+1,2*j  )) >= 1.e-15){
                                int dof_i = 3*numNodesCoarse + 2 * connec(i) +1;
                                int dof_j = 3*numNodesCoarse + 2 * connec(j);
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&localMV.first(2*i+1,2*j  ),ADD_VALUES);
                            };

                            if (fabs(localMV.first(2*i  ,2*j+1)) >= 1.e-15){
                                int dof_i = 3*numNodesCoarse + 2 * connec(i);
                                int dof_j = 3*numNodesCoarse + 2 * connec(j) +1;
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&localMV.first(2*i  ,2*j+1),ADD_VALUES);
                            };

                            if (fabs(localMV.first(2*i+1,2*j+1)) >= 1.e-15){
                                int dof_i = 3*numNodesCoarse + 2 * connec(i) +1;
                                int dof_j = 3*numNodesCoarse + 2 * connec(j) +1;
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&localMV.first(2*i+1,2*j+1),ADD_VALUES);
                            };

                            //Matrix Q and Qt
                            if (fabs(localMV.first(2*i  ,12+j)) >= 1.e-15){
                                int dof_i = 3 * numNodesCoarse + 2 * connec(i);
                                int dof_j = 3 * numNodesCoarse + 2 * numNodesFine + connec(j);
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&localMV.first(2*i  ,12+j),ADD_VALUES);
                            };

                            if (fabs(localMV.first(12+j,2*i  )) >= 1.e-15){
                                int dof_i = 3 * numNodesCoarse + 2 * connec(i);
                                int dof_j = 3 * numNodesCoarse + 2 * numNodesFine + connec(j);
                                ierr = MatSetValues(A,1,&dof_j,1,&dof_i,&localMV.first(12+j,2*i  ),ADD_VALUES);
                            };

                            if (fabs(localMV.first(2*i+1,12+j)) >= 1.e-15){
                                int dof_i = 3 * numNodesCoarse + 2 * connec(i) + 1;
                                int dof_j = 3 * numNodesCoarse + 2 * numNodesFine + connec(j);
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&localMV.first(2*i+1,12+j),ADD_VALUES);
                            };
                            
                            if (fabs(localMV.first(12+j,2*i+1)) >= 1.e-15){
                                int dof_i = 3 * numNodesCoarse + 2 * connec(i) + 1;
                                int dof_j = 3 * numNodesCoarse + 2 * numNodesFine + connec(j);
                                ierr = MatSetValues(A,1,&dof_j,1,&dof_i,&localMV.first(12+j,2*i+1),ADD_VALUES);
                            };
                            
                            if (fabs(localMV.first(12+i,12+j)) >= 1.e-15){
                                int dof_i = 3 * numNodesCoarse + 2 * numNodesFine + connec(i);
                                int dof_j = 3 * numNodesCoarse + 2 * numNodesFine + connec(j);
                                ierr = MatSetValues(A,1,&dof_j,1,&dof_i,&localMV.first(12+i,12+j),ADD_VALUES);
                            };
                        };     
                        ///Rhs vector
                        if (fabs(localMV.second(2*i  )) >= 1.e-15){
                            int dof_i = 3 * numNodesCoarse + 2 * connec(i);
                            ierr = VecSetValues(b,1,&dof_i,&localMV.second(2*i  ),ADD_VALUES);
                        };

                        if (fabs(localMV.second(2*i+1)) >= 1.e-15){
                            int dof_i = 3 * numNodesCoarse + 2 * connec(i) + 1;
                            ierr = VecSetValues(b,1,&dof_i,&localMV.second(2*i+1),ADD_VALUES);
                        };

                        if (fabs(localMV.second(12+i)) >= 1.e-15){
                            int dof_i = 3 * numNodesCoarse + 2 * numNodesFine + connec(i);
                            ierr = VecSetValues(b,1,&dof_i,&localMV.second(12+i),ADD_VALUES);
                        };
                    }; 
                };                
            };
              

            //Lagrange Multipliers
            for (int l=0; l< numElemGlueZoneFine; l++){
                
                int jel = elementsGlueZoneFine_[l];
                
                if (domDecompFine.first[jel] == rank) {
                    
                    typename Elements::LocalMatrix Ajac, AjacAnt, AStab, ArlequinM, ArlequinM2;
                    typename Elements::LocalVector Rhs, RhsStab, RhsArlequin,rhsLagMult;
                    typename Elements::Connectivity connec, connecC, connecL;
                    
                    connec = elementsFine_[jel] -> getConnectivity();
                    connecL = glueZoneFine_[l] -> getConnectivity();

                    // FINE MESH
                    
                    //Computes element matrix
                    elementsFine_[jel] -> getLagrangeMultipliersSameMesh();
                    
                    //Gets element matrice and rhs vectors
                    // -L1
                    Ajac = - elementsFine_[jel] -> getLagrMultMatrix();
                    // +L1 * Lambda
                    // rhsLagMult = - elementsFine_[jel] -> getRhsLagrangeMultipliers(-Ajac);
                    rhsLagMult = - elementsFine_[jel] -> getLagrMultVector();

                    //PSPG and SUPG stabilizations
                    AStab = elementsFine_[jel] -> getJacNRMatrix();
                    RhsStab = elementsFine_[jel] -> getRhsVector();

                    //Arlequin Stabilization
                    ArlequinM = elementsFine_[jel] -> getArlequinStabilizationMatrix();
                    RhsArlequin = elementsFine_[jel] -> getArlequinStabilizationVector();

                    ArlequinM2 = elementsFine_[jel] -> getArlequinStabilizationMatrix2();

                    // -L1t
                    Ajac = trans(Ajac);

                    // L1t * u1
                    Rhs = -elementsFine_[jel] -> getRhsVelocities(-Ajac);

                    std::pair<Elements::LocalVector, Elements::LocalMatrix>  lagMult;
   
                    lagMult.first.clear();
                    lagMult.second.clear();
                         
                    lagMult = elementsFine_[jel] -> getBoundaryConditionsVelocity(Rhs,Ajac);

                    Ajac = lagMult.second;
                    Rhs = lagMult.first;


                    lagMult.first.clear();
                    lagMult.second.clear();
                         
                    lagMult = elementsFine_[jel] -> getBoundaryConditionsLagMult(rhsLagMult,Ajac);

                    Ajac = lagMult.second;
                    rhsLagMult = lagMult.first;


                    double integ = alpha_f * gamma * dTime;
                    //Disperse local contributions into the global matrix
                    for (int i=0; i<6; i++){
                        for (int j=0; j<6; j++){
                            
                            //COUPLING OPERATOR
                            if (fabs(Ajac(2*i  ,2*j  )) >= 1.e-15){
                                int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                                int d_j = 3*numNodesCoarse + 2*connec(j);
                                double value = Ajac(2*i  ,2*j  )*integ;
                                ierr = MatSetValues(A,1,&d_i,1,&d_j,&value,ADD_VALUES);
                                ierr = MatSetValues(A,1,&d_j,1,&d_i,&Ajac(2*i  ,2*j  ),ADD_VALUES);
                            };
                            if (fabs(Ajac(2*i+1,2*j  )) >= 1.e-15){
                                int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i) + 1;
                                int d_j = 3*numNodesCoarse + 2*connec(j);
                                double value = Ajac(2*i+1,2*j  ) * integ;
                                ierr = MatSetValues(A,1,&d_i,1,&d_j,&value,ADD_VALUES);
                                ierr = MatSetValues(A,1,&d_j,1,&d_i,&Ajac(2*i+1,2*j  ),ADD_VALUES);
                            };
                            if (fabs(Ajac(2*i+1,2*j+1)) >= 1.e-15){
                                int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i) + 1;
                                int d_j = 3*numNodesCoarse + 2*connec(j) + 1;
                                double value = Ajac(2*i+1,2*j+1)*integ;
                                ierr = MatSetValues(A,1,&d_i,1,&d_j,&value,ADD_VALUES);
                                ierr = MatSetValues(A,1,&d_j,1,&d_i,&Ajac(2*i+1,2*j+1),ADD_VALUES);
                            };
                            if (fabs(Ajac(2*i  ,2*j+1)) >= 1.e-15){
                                int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                                int d_j = 3 * numNodesCoarse + 2*connec(j) + 1;
                                double value = Ajac(2*i  ,2*j+1)*integ;
                                ierr = MatSetValues(A,1,&d_i,1,&d_j,&value,ADD_VALUES);
                                ierr = MatSetValues(A,1,&d_j,1,&d_i,&Ajac(2*i  ,2*j+1),ADD_VALUES);
                            };

                            //SUPG STABILIZATION
                            if (fabs(AStab(2*i  ,2*j  )) >= 1.e-15){
                                int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                                int d_j = 3*numNodesCoarse + 2*connec(j);
                                ierr = MatSetValues(A,1,&d_j,1,&d_i,&AStab(2*i  ,2*j  ),ADD_VALUES);
                            };
                            if (fabs(AStab(2*i+1,2*j+1)) >= 1.e-15){
                                int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i) + 1;
                                int d_j = 3*numNodesCoarse + 2*connec(j) + 1;
                                ierr = MatSetValues(A,1,&d_j,1,&d_i,&AStab(2*i+1,2*j+1),ADD_VALUES);
                            };

                            //PSPG STABILIZATION
                            if (fabs(AStab(2*i  ,12+j)) >= 1.e-15){
                                int dof_i = 3*numNodesCoarse + 2*numNodesFine + connec(i);
                                int dof_j = 3*numNodesCoarse + 3*numNodesFine + connecL(j);
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&AStab(2*i  ,12+j),ADD_VALUES);
                            };
                            if (fabs(AStab(2*i+1,12+j)) >= 1.e-15){
                                int dof_i = 3*numNodesCoarse + 2*numNodesFine + connec(i);
                                int dof_j = 3*numNodesCoarse + 3*numNodesFine + connecL(j)+1;
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&AStab(2*i+1,12+j),ADD_VALUES);
                            };

                            //ARLEQUIN STABILIZATION
                            // if (fabs(ArlequinM(2*i  ,2*j+1)) >= 1.e-15){
                            //     int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                            //     int d_j = 3*numNodesCoarse + 2*connec(j) + 1;
                            //     ierr = MatSetValues(A,1,&d_i,1,&d_j,
                            //                         &ArlequinM(2*i  ,2*j+1),ADD_VALUES);
                            // };
                            if (fabs(ArlequinM2(2*i  ,2*j  )) >= 1.e-15){
                                int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                                int d_j = 3*numNodesCoarse + 2*connec(j);
                                ierr = MatSetValues(A,1,&d_i,1,&d_j,&ArlequinM2(2*i  ,2*j  ),ADD_VALUES);
                            };
                            // if (fabs(ArlequinM(2*i+1,2*j  )) >= 1.e-15){
                            //     int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i) + 1;
                            //     int d_j = 3*numNodesCoarse + 2*connec(j);
                            //     ierr = MatSetValues(A,1,&d_i,1,&d_j,
                            //                         &ArlequinM(2*i+1,2*j  ),ADD_VALUES);
                            // };
                            if (fabs(ArlequinM2(2*i+1,2*j+1)) >= 1.e-15){
                                int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i) + 1;
                                int d_j = 3*numNodesCoarse + 2*connec(j) + 1;
                                ierr = MatSetValues(A,1,&d_i,1,&d_j,&ArlequinM2(2*i+1,2*j+1),ADD_VALUES);
                            };

                            if (fabs(ArlequinM2(12+i,2*j  )) >= 1.e-15){
                                int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                                int dof_j = 3*numNodesCoarse + 2*numNodesFine + connec(j);
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinM2(12+i,2*j  ),ADD_VALUES);
                            };
                            if (fabs(ArlequinM2(12+i,2*j+1)) >= 1.e-15){
                                int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i) + 1;
                                int dof_j = 3*numNodesCoarse + 2*numNodesFine + connec(j);
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinM2(12+i,2*j+1),ADD_VALUES);
                            };
                            if (fabs(ArlequinM(2*i  ,2*j  )) >= 1.e-15){
                                int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                                int dof_j = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(j);
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinM(2*i  ,2*j  ),ADD_VALUES);
                                dof_i++; dof_j++;
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinM(2*i+1,2*j+1),ADD_VALUES);
                            };                            
                        };

                        //RHS VECTOR
                        //COUPLING OPERATOR
                        int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                        ierr = VecSetValues(b,1,&dof_i,&Rhs(2*i  ),ADD_VALUES);

                        dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i) + 1;
                        ierr = VecSetValues(b,1,&dof_i,&Rhs(2*i+1),ADD_VALUES);

                        dof_i = 3*numNodesCoarse + 2*connec(i);
                        ierr = VecSetValues(b,1,&dof_i,&rhsLagMult(2*i  ),ADD_VALUES);

                        dof_i = 3*numNodesCoarse + 2*connec(i) + 1;
                        ierr = VecSetValues(b,1,&dof_i,&rhsLagMult(2*i+1),ADD_VALUES);

                        //SUPG STABILIZATION
                        dof_i = 3*numNodesCoarse + 2*connec(i);                       
                        ierr = VecSetValues(b,1,&dof_i,&RhsStab(2*i  ),ADD_VALUES);

                        dof_i = 3*numNodesCoarse + 2*connec(i) + 1;
                        ierr = VecSetValues(b,1,&dof_i,&RhsStab(2*i+1),ADD_VALUES);

                        //PSPG STABILIZATION
                        dof_i = 3*numNodesCoarse + 2*numNodesFine + connec(i);
                        ierr = VecSetValues(b,1,&dof_i,&RhsStab(12+i),ADD_VALUES);

                        // ///ARLEQUIN STABILIZATION
                        dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                        ierr = VecSetValues(b,1,&dof_i,&RhsArlequin(2*i  ),ADD_VALUES);
                        dof_i++;
                        ierr = VecSetValues(b,1,&dof_i,&RhsArlequin(2*i+1),ADD_VALUES);
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
                        double pspg = 0;//elementsCoarse_[iElemCoarse] -> getPSPG();
                        ublas::bounded_vector<double, 6> press_, velX_, velY_;

                        connecC = elementsCoarse_[iElemCoarse] -> 
                            getConnectivity();

                        for (int k = 0; k < 6; k++){
                            press_(k) = nodesCoarse_[connecC(k)] -> getPressure();
                            velX_(k) = nodesCoarse_[connecC(k)] -> getVelocity(0);
                            velY_(k) = nodesCoarse_[connecC(k)] -> getVelocity(1);
                        }
                        
                        elementsFine_[jel] -> 
                            getLagrangeMultipliersDifferentMesh(iElemCoarse,pspg,press_,velX_,velY_);
                        
                        //Computes element matrix
                        // L0
                        Ajac = elementsFine_[jel] -> getLagrMultMatrix();
                        // L0t
                        Ajac = trans(Ajac);
                        AStab = elementsFine_[jel] -> getJacNRMatrix();
                        AStab = trans(AStab);
                        RhsStab = elementsFine_[jel] -> getRhsVector();
                        ArlequinM = elementsFine_[jel] -> getArlequinStabilizationMatrix();
                        ArlequinM = trans(ArlequinM);
                        RhsArlequin = elementsFine_[jel] -> getArlequinStabilizationVector();
                        
                        ArlequinM2 = elementsFine_[jel] -> getArlequinStabilizationMatrix2();
                        ArlequinM2 = trans(ArlequinM2);

                        // -L0 * lambda
                        rhsLagMult =  elementsFine_[jel] -> getRhsLagrangeMultipliers(trans(Ajac));

                        

                        // AStab.clear();
                        // // RhsStab.clear();
                        // ArlequinM.clear();
                        // RhsArlequin.clear();


                        
                        // U_.clear();
                        // for (int i = 0; i < 6; i++){
                        //     U_(2*i  ) = nodesFine_[connec(i)] -> getLagrangeMultiplier(0);
                        //     U_(2*i+1) = nodesFine_[connec(i)] -> getLagrangeMultiplier(1);
                        // };
                        
                        // noalias(rhsLagMult) = -prod(trans(Ajac),U_);
                        // noalias(lagStab) = -prod(trans(AjacAnt),U_);
                        


 
                    Rhs = elementsCoarse_[iElemCoarse] -> getRhsVelocities(Ajac);

                    std::pair<Elements::LocalVector, Elements::LocalMatrix>  lagMult;
                    


                    // lagMult.first.clear();
                    // lagMult.second.clear();
                         
                    // lagMult = elementsCoarse_[iElemCoarse] -> getBoundaryConditionsVelocity(rhsLagMult,(Ajac));

                    // Ajac = lagMult.second;
                    // rhsLagMult = lagMult.first;


                    // lagMult.first.clear();
                    // lagMult.second.clear();
                         
                    // lagMult = elementsFine_[jel] -> getBoundaryConditionsLagMult(Rhs,(Ajac));

                    // Ajac = (lagMult.second);
                    // Rhs = lagMult.first;


                    // if (jel == 94){
                    //     std::cout << "MATRIZ " << std::endl;
                    //     for (int i=0; i<18; i++){
                    //         for (int j=0; j<18; j++){
                    //             std::cout << Ajac(i,j) << " ";
                    //         }
                    //         std::cout << std::endl;
                    //     }
                    //     std::cout << "VETOR " << std::endl;
                    //     for (int i=0; i<18; i++){
                    //         std::cout << Rhs(i) << " " << rhsLagMult(i) << std::endl;
                    //     }
                    // }

                    rhsLagMult *= 1000.e0;
                    //Ajac *= 1000.;

                    for (int i = 0; i < 12; i++){
                        RhsStab(i) *= 1000.e0;
                        for (int j = 0; j < 12; j++){
                            AStab(i,j)*=1000.e0;
                            Ajac(i,j)*=1000.e0;
                            ArlequinM2(i,j)*=1000.e0;
                        }
                    }
                        // std::pair<Elements::LocalVector, 
                        //           Elements::LocalMatrix>  lagMult;
                        
                        // lagMult.first.clear();
                        // lagMult.second.clear();
                         
                        // lagMult = elementsCoarse_[iElemCoarse] -> 
                        //     getRhsVectorAndBoundaryConditions(Ajac);

                        // noalias(rhsLagMult) = -prod((lagMult.second),U_);
                    // for (int i=0; i<6; i++){
                    //     if(lagMult.second(12+i,12+i) > 0){
                    //         rhsLagMult(2*i  ) = 0.;
                    //         rhsLagMult(2*i+1) = 0.;
                    //         double one = 1.;
                    //         int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                    //         ierr = MatSetValues(A,1,&d_i,1,&d_i,&one,INSERT_VALUES);
                    //         d_i++;
                    //         ierr = MatSetValues(A,1,&d_i,1,&d_i,&one,INSERT_VALUES);
                    //     }
                    // }
                        double integ = alpha_f * gamma * dTime;
                        //Disperse local contribution into the global matrix
                        for (int i=0; i<6; i++){
                            for (int j=0; j<6; j++){
                                //COUPLING OPERATOR
                                if (fabs(Ajac(2*i  ,2*j  )) >= 1.e-15){
                                    int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                                    int dof_j = 2*connecC(j);
                                    double value = Ajac(2*i  ,2*j  ) * integ;
                                    ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&value,ADD_VALUES);
                                    ierr = MatSetValues(A,1,&dof_j,1,&dof_i,&Ajac(2*i  ,2*j  ),ADD_VALUES);
                                }
                                if (fabs(Ajac(2*i+1,2*j  )) >= 1.e-15){
                                    int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i) + 1;
                                    int dof_j = 2*connecC(j);
                                    double value = Ajac(2*i+1,2*j  )*integ;
                                    ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&value,ADD_VALUES);
                                    ierr = MatSetValues(A,1,&dof_j,1,&dof_i,&Ajac(2*i+1,2*j  ),ADD_VALUES);
                                };
                                if (fabs(Ajac(2*i+1,2*j+1)) >= 1.e-15){
                                    int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i) + 1;
                                    int dof_j = 2*connecC(j) + 1;
                                    double value = Ajac(2*i+1,2*j+1)*integ;
                                    ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&value,ADD_VALUES);
                                    ierr = MatSetValues(A,1,&dof_j,1,&dof_i,&Ajac(2*i+1,2*j+1),ADD_VALUES);
                                };
                                if (fabs(Ajac(2*i  ,2*j+1)) >= 1.e-15){
                                    int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                                    int dof_j = 2*connecC(j) + 1;
                                    double value = Ajac(2*i  ,2*j+1)*integ;
                                    ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&value,ADD_VALUES);
                                    ierr = MatSetValues(A,1,&dof_j,1,&dof_i,&Ajac(2*i  ,2*j+1),ADD_VALUES);
                                };
                                if (fabs(Ajac(12+i,12+j)) >= 1.e-15){
                                    int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                                    int dof_j = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(j);
                                    // ierr = MatSetValues(A,1,&dof_i,1,&dof_j,
                                    //                 &lagMult.second(12+i,12+j),ADD_VALUES);
                                    dof_i++; dof_j++;
                                    ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&Ajac(12+i,12+j),ADD_VALUES);
                                };
                         

                                //SUPG STABILIZATION
                                if (fabs(AStab(2*i  ,2*j  )) >= 1.e-15){
                                    int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                                    int dof_j = 2*connecC(j);
                                    ierr = MatSetValues(A,1,&dof_j,1,&dof_i,&AStab(2*i  ,2*j  ),ADD_VALUES);
                                };
                                if (fabs(AStab(2*i+1,2*j+1)) >= 1.e-15){
                                    int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i) + 1;
                                    int dof_j = 2*connecC(j) + 1;
                                    ierr = MatSetValues(A,1,&dof_j,1,&dof_i,&AStab(2*i+1,2*j+1),ADD_VALUES);
                                };

                                //PSPG STABILIZATION
                                if (fabs(AStab(2*i  ,12+j)) >= 1.e-15){
                                    int dof_i = 2*numNodesCoarse + connecC(i);
                                    int dof_j = 3*numNodesCoarse + 2*numNodesFine + connecL(j);
                                    ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&AStab(2*i  ,12+j),ADD_VALUES);
                                };
                                if (fabs(AStab(2*i+1,12+j)) >= 1.e-15){
                                    int dof_i = 2*numNodesCoarse + connecC(i);
                                    int dof_j = 3*numNodesCoarse + 2*numNodesFine + connecL(j)+1;
                                    ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&AStab(2*i+1,12+j),ADD_VALUES);
                                };

                                ///ARLEQUIN STABILIZATION
                                // if (fabs(ArlequinM(2*i  ,2*j+1)) >= 1.e-15){
                                //     int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                                //     int d_j = 2*connecC(j) + 1;
                                //     ierr = MatSetValues(A,1,&d_i,1,&d_j,
                                //                         &ArlequinM(2*i  ,2*j+1),ADD_VALUES);
                                // };
                                if (fabs(ArlequinM2(2*i  ,2*j  )) >= 1.e-15){
                                    int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                                    int d_j = 2*connecC(j);
                                    ierr = MatSetValues(A,1,&d_i,1,&d_j,&ArlequinM2(2*i  ,2*j  ),ADD_VALUES);
                                };
                                // if (fabs(ArlequinM(2*i+1,2*j  )) >= 1.e-15){
                                //     int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i) + 1;
                                //     int d_j = 2*connecC(j);
                                //     ierr = MatSetValues(A,1,&d_i,1,&d_j,
                                //                         &ArlequinM(2*i+1,2*j  ),ADD_VALUES);
                                // };
                                if (fabs(ArlequinM2(2*i+1,2*j+1)) >= 1.e-15){
                                    int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i) + 1;
                                    int d_j = 2*connecC(j) + 1;
                                    ierr = MatSetValues(A,1,&d_i,1,&d_j,&ArlequinM2(2*i+1,2*j+1),ADD_VALUES);
                                };

                                if (fabs(ArlequinM2(12+i,2*j  )) >= 1.e-15){
                                    int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                                    int dof_j = 2*numNodesCoarse + connecC(j);
                                    ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinM2(12+i,2*j  ),ADD_VALUES);
                                };
                                if (fabs(ArlequinM2(12+i,2*j+1)) >= 1.e-15){
                                    int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i) + 1;
                                    int dof_j = 2*numNodesCoarse + connecC(j);
                                    ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinM2(12+i,2*j+1),ADD_VALUES);
                                };

                                if (fabs(ArlequinM(2*i  ,2*j  )) >= 1.e-15){
                                    int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                                    int dof_j = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(j);
                                    ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinM(2*i  ,2*j  ),ADD_VALUES);
                                    dof_i++; dof_j++;
                                    ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinM(2*i+1,2*j+1),ADD_VALUES);
                                };
                            };
                            //RHS VECTOR
                            //COUPLING OPERATOR
                            int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                            ierr = VecSetValues(b,1,&d_i,&Rhs(2*i  ),ADD_VALUES);

                            d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i) + 1;
                            ierr = VecSetValues(b,1,&d_i,&Rhs(2*i+1),ADD_VALUES);

                            int dof_i = 2*connecC(i);
                            ierr = VecSetValues(b,1,&dof_i,&rhsLagMult(2*i  ),ADD_VALUES);

                            dof_i = 2*connecC(i) + 1;
                            ierr = VecSetValues(b,1,&dof_i,&rhsLagMult(2*i+1),ADD_VALUES);

                            //SUPG STABILIZATION
                            dof_i = 2*connecC(i);
                            ierr = VecSetValues(b,1,&dof_i,&RhsStab(2*i  ),ADD_VALUES);

                            dof_i = 2*connecC(i) + 1;
                            ierr = VecSetValues(b,1,&dof_i,&RhsStab(2*i+1),ADD_VALUES);

                            //PSPG STABILIZATION
                            dof_i = 2*numNodesCoarse + connecC(i);
                            ierr = VecSetValues(b,1,&dof_i,&RhsStab(12+i),ADD_VALUES);

                            //ARLEQUIN STABILIZATION
                            dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                            ierr = VecSetValues(b,1,&dof_i,&RhsArlequin(2*i  ),ADD_VALUES);
                            dof_i++;
                            ierr = VecSetValues(b,1,&dof_i,&RhsArlequin(2*i+1),ADD_VALUES);
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


            // MatScale(A,1000.);
            // VecScale(b,1000.);


            // Vec AuxRow, AuxCol;
            // ierr = VecDuplicate(b, &AuxRow); CHKERRQ(ierr);
            // ierr = VecDuplicate(b, &AuxCol); CHKERRQ(ierr);
            // Vec v1, v2;
            // MatCreateVecs(A,&v1,&v2);
            // PetscInt idxm[sysSize];
            // for (int i = 0; i < sysSize; ++i){
            //     idxm[i] = i;
            // }

            // for (int i = 0; i < numNodesCoarse; ++i){
            //     if (nodesCoarse_[i] -> getWeightFunction() < 1){
                    
            //         int lin = 2*i;
            //         MatGetColumnVector(A,AuxCol,lin);
            //         VecScale(AuxCol, 1000);
            //         MatGetRow(A,AuxRow,NULL,NULL,NULL);
            //         VecScale(AuxRow, 1000);
            //         MatZeroRowsColumns(A,1,&lin,0.,&v1,&v2);
            //         MatSetValues(A,1,&lin,sysSize,idxm,AuxRow,ADD_VALUES);
            //         MatSetValues(A,sysSize,idxm,1,&lin,AuxCol,ADD_VALUES);
            //     }
            // }

            // ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
            // ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

            // VecDestroy(&v1);
            // VecDestroy(&v2);
            // VecDestroy(&AuxRow);
            // VecDestroy(&AuxCol);

            
 // PetscViewer    viewer;

 // PetscViewerDrawOpen(PETSC_COMM_WORLD,NULL,NULL,0,0,300,300,&viewer);
 // PetscObjectSetName((PetscObject)viewer,"Line graph Plot");
 //  PetscViewerPushFormat(viewer,PETSC_VIEWER_DRAW_LG);

            //ierr = MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
            //ierr = VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
            
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
            // ierr = KSPSetType(ksp,KSPGMRES); CHKERRQ(ierr);
            
            // ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
            // // ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);



            
#if defined(PETSC_HAVE_MUMPS)
            ierr = KSPSetType(ksp,KSPPREONLY);
            ierr = KSPGetPC(ksp,&pc);
            ierr = PCSetType(pc, PCLU);
            
            ierr = PCFactorSetMatSolverType(pc,MATSOLVERMUMPS);
            PCFactorSetUpMatSolverType(pc);
            PCFactorGetMatrix(pc,&F);
            
            PetscInt ival,icntl;
            icntl = 14; ival = 80;
            MatMumpsSetIcntl(F,icntl,ival);
            icntl = 28; ival = 2;
            MatMumpsSetIcntl(F,icntl,ival);
            icntl = 29; ival = 2;
            MatMumpsSetIcntl(F,icntl,ival);
            icntl = 16; ival = 0;
            MatMumpsSetIcntl(F,icntl,ival);
            // icntl = 4; ival = 3;
            // MatMumpsSetIcntl(F,icntl,ival);
            // icntl = 11; ival = 1;
            // MatMumpsSetIcntl(F,11,1);

            //MatMumpsSetIcntl(F,21,0);

            
#endif
            ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
            ierr = KSPSetUp(ksp);
            
#if defined(PETSC_HAVE_MUMPS)
            PetscInt  info1,info2,icntl14;
            
            MatMumpsGetInfo(F,1,&info1);
            MatMumpsGetInfo(F,2,&info2);

            // PetscReal info5,info6,info7,info8,info9,info10,info11;
            // MatMumpsGetRinfog(F,5,&info5);
            // MatMumpsGetRinfog(F,6,&info6);
            // MatMumpsGetRinfog(F,7,&info7);
            // MatMumpsGetRinfog(F,8,&info8);
            // MatMumpsGetRinfog(F,9,&info9);
            // MatMumpsGetRinfog(F,10,&info10);
            // MatMumpsGetRinfog(F,11,&info11);

            // PetscInt info21,info32;
            // MatMumpsGetIcntl(F,32,&info32);
            // MatMumpsGetIcntl(F,21,&info21);

            
            // if(rank==0) std::cout << "ICNTL = " << info21 << " " << info32 << std::endl;
          
            // if(rank==0) std::cout << "INFOG = " << info5 << " " << info6 << " " << info7 << " " << info8 << " " << info9 << " " << info10 << " " << info11 << std::endl;
            MatMumpsGetIcntl(F,14,&icntl14);    
            if((rank == 0) && (info1 != 0)) std::cout << " INFO(1) = " << info1
                                                      << " " << info2 << " " 
                                                      << icntl14 << std::endl;
#endif
            
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
            double normU = 0.;
            double normP = 0.;
            double normL = 0.;
            double normT = 0.;
            double p_;
            Ione = 1;

            


            for (int i = 0; i < numNodesCoarse; ++i){
                Ii = 2 * i;
                ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                if (nodesCoarse_[i] -> getDistFunction() > -1.2) val *= 1000.e0;
                u_[0] = val;
                normU += val*val;
                nodesCoarse_[i] -> incrementAcceleration(0,u_[0]);
                nodesCoarse_[i] -> incrementVelocity(0,u_[0]*gamma*dTime);
                
                Ii = 2 * i + 1;
                ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                if (nodesCoarse_[i] -> getDistFunction() > -1.2) val *= 1000.e0;
                u_[1] = val;
                normU += val*val;
                nodesCoarse_[i] -> incrementAcceleration(1,u_[1]);
                nodesCoarse_[i] -> incrementVelocity(1,u_[1]*gamma*dTime);
            };
            
            for (int i = 0; i<numNodesCoarse; i++){
                Ii = 2 * numNodesCoarse + i;
                ierr = VecGetValues(All,Ione,&Ii,&val);CHKERRQ(ierr);
                p_ = val;
                normP += val*val;
                nodesCoarse_[i] -> incrementPressure(p_);
            };
            
            for (int i = 0; i < numNodesFine; ++i){
                Ii = 3 * numNodesCoarse + 2 * i;
                ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                u_[0] = val;
                normU += val*val;
                nodesFine_[i] -> incrementAcceleration(0,u_[0]);
                nodesFine_[i] -> incrementVelocity(0,u_[0]*gamma*dTime);
                
                Ii = 3 * numNodesCoarse + 2 * i + 1;
                ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                u_[1] = val;
                normU += val*val;
                nodesFine_[i] -> incrementAcceleration(1,u_[1]);
                nodesFine_[i] -> incrementVelocity(1,u_[1]*gamma*dTime);
            };
            
            for (int i = 0; i<numNodesFine; i++){
                Ii = 3 * numNodesCoarse + 2 * numNodesFine + i;
                ierr = VecGetValues(All,Ione,&Ii,&val);CHKERRQ(ierr);
                p_ = val;
                normP += val*val;
                nodesFine_[i] -> incrementPressure(p_);
            };
            
            for (int i = 0; i < numNodesGlueZoneFine; ++i){
                Ii = 3 * numNodesCoarse + 3 * numNodesFine + 2 * i;
                ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                u_[0] = val;
                nodesFine_[nodesGlueZoneFine_[i]] -> 
                    incrementLagrangeMultiplier(0,u_[0]);
                normL += val*val;

                Ii = 3 * numNodesCoarse + 3 * numNodesFine + 2 * i + 1;
                ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                u_[1] = val;
                nodesFine_[nodesGlueZoneFine_[i]] -> 
                    incrementLagrangeMultiplier(1,u_[1]);
                normL += val*val;
                // std::cout << "LAG M " << u_[0] << " " << u_[1] << std::endl;
            };
            
            //Computes the solution vector norm
            ierr = VecNorm(u,NORM_2,&val);CHKERRQ(ierr);
            
            boost::posix_time::ptime t2 =
                boost::posix_time::microsec_clock::local_time();
            
            if(rank == 0){
                boost::posix_time::time_duration diff = t2 - t1;
                
                std::cout<<"Iteration = " << inewton << " (" << iterations <<  
                    ")  Du Norm = " << std::scientific << sqrt(normU) 
                         << " " << sqrt(normP) 
                         << " " << sqrt(normL) 
                         << " " << val << 
                    "  Time (s) = " << std::fixed << 
                    diff.total_milliseconds()/1000. << std::endl;
            };
            
            ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
            ierr = VecDestroy(&b); CHKERRQ(ierr);
            ierr = VecDestroy(&u); CHKERRQ(ierr);
            ierr = VecDestroy(&All); CHKERRQ(ierr);
            ierr = MatDestroy(&A); CHKERRQ(ierr);
           // ierr = MatDestroy(&F); CHKERRQ(ierr);
            //ierr = MatDestroy(&C); CHKERRQ(ierr);
            
            if(val <= tolerance){
                break;            
            };          
            
        };
        

        if (rank == 0){
            double normUUprev, normU;
            normUUprev = 0.;
            normU = 0.;
            for (int i = 0; i < numNodesCoarse; ++i){
                double u, v, uPr, vPr, weight;
                weight = nodesCoarse_[i] -> getWeightFunction();
                u = nodesCoarse_[i] -> getVelocity(0) * weight;
                v = nodesCoarse_[i] -> getVelocity(1) * weight;
                uPr = nodesCoarse_[i] -> getPreviousVelocity(0) * weight;
                vPr = nodesCoarse_[i] -> getPreviousVelocity(1) * weight;

                normUUprev += (u-uPr) * (u-uPr) + (v-vPr) * (v-vPr);
                normU += u*u + v*v;
            }
            for (int i = 0; i < numNodesFine; ++i){
                double u, v, uPr, vPr, weight;
                weight = nodesFine_[i] -> getWeightFunction();
                u = nodesFine_[i] -> getVelocity(0) * weight;
                v = nodesFine_[i] -> getVelocity(1) * weight;
                uPr = nodesFine_[i] -> getPreviousVelocity(0) * weight;
                vPr = nodesFine_[i] -> getPreviousVelocity(1) * weight;

                normUUprev += (u-uPr) * (u-uPr) + (v-vPr) * (v-vPr);
                normU += u*u + v*v;
            }
            std::cout << "NORM U  " << std::scientific <<  sqrt(normUUprev / normU) << std::endl;
        }



        //Compute real velocity
        QuadShapeFunction<2>                       shapeQuad;
        typename QuadShapeFunction<2>::Values      phi_;
        
        for (int i = 0; i<numNodesFine; i++){
            nodesFine_[i] -> setVelocityArlequin(0,nodesFine_[i] -> 
                                                 getVelocity(0));
            nodesFine_[i] -> setVelocityArlequin(1,nodesFine_[i] -> 
                                                 getVelocity(1));
            nodesFine_[i] -> setPressureArlequin(nodesFine_[i] ->getPressure());
        };
        
        for (int i = 0; i<numNodesGlueZoneFine; i++){
            typename Nodes::VecLocD xsi;
            typename Quadrature::NodalValuesQuad u_coarse, v_coarse, p_coarse;
            typename Elements::Connectivity connecCoarse;
            
            double u = 0.;
            double v = 0.;
            double p = 0.;
            
            int elCoarse = nodesFine_[nodesGlueZoneFine_[i]] -> 
                getNodalElemCorrespondence();
            xsi = nodesFine_[nodesGlueZoneFine_[i]] -> 
                getNodalXsiCorrespondence();
            
            connecCoarse = elementsCoarse_[elCoarse] -> getConnectivity();
            
            for (int j=0; j<6; j++){
                u_coarse(j) = nodesCoarse_[connecCoarse(j)] -> getVelocity(0);
                v_coarse(j) = nodesCoarse_[connecCoarse(j)] -> getVelocity(1);
                p_coarse(j) = nodesCoarse_[connecCoarse(j)] -> getPressure();
            };
            
            shapeQuad.evaluate(xsi,phi_);
            
            for (int j=0; j<6; j++){
                u += u_coarse(j) * phi_(j);
                v += v_coarse(j) * phi_(j);
                p += p_coarse(j) * phi_(j);
            };
            
            double wFunc = nodesFine_[nodesGlueZoneFine_[i]] -> 
                getWeightFunction();
            
            double u_int = nodesFine_[nodesGlueZoneFine_[i]] -> 
                getVelocity(0) * wFunc + u * (1. - wFunc);
            double v_int = nodesFine_[nodesGlueZoneFine_[i]] -> 
                getVelocity(1) * wFunc + v * (1. - wFunc);
            double p_int = nodesFine_[nodesGlueZoneFine_[i]] -> 
                getPressure() * wFunc + p * (1. - wFunc);
            
            nodesFine_[nodesGlueZoneFine_[i]] -> setVelocityArlequin(0,u_int);
            nodesFine_[nodesGlueZoneFine_[i]] -> setVelocityArlequin(1,v_int);
            nodesFine_[nodesGlueZoneFine_[i]] -> setPressureArlequin(p_int);
            
        };
        
        for (int i=0; i<numNodesCoarse; i++){
            nodesCoarse_[i] -> 
                setVelocityArlequin(0,nodesCoarse_[i] -> getVelocity(0));
            nodesCoarse_[i] -> 
                setVelocityArlequin(1,nodesCoarse_[i] -> getVelocity(1));
            nodesCoarse_[i] -> 
                setPressureArlequin(nodesCoarse_[i] -> getPressure());
        };
        
        // Compute and print drag and lift coefficients
        if (fineModel.getComputeDragAndLift()){
            dragAndLiftCoefficients(dragLift);
        };

        if (rank == 0) {
            if (coarseModel.printVorticity){
                for (int i = 0; i < numNodesCoarse; i++){
                    nodesCoarse_[i] -> clearVorticity();
                };
                for (int i = 0; i < numNodesFine; i++){
                    nodesFine_[i] -> clearVorticity();
                };
                for (int jel = 0; jel < numElemCoarse; jel++){
                    elementsCoarse_[jel] -> computeVorticity();
                };
                for (int jel = 0; jel < numElemFine; jel++){
                    elementsFine_[jel] -> computeVorticity();
                };
            };

            //Printing results
            printResults(iTimeStep);
        };        
 

    };
        
    return 0;

};



//------------------------------------------------------------------------------
//----------------COMPUTE ARLEQUIN COUPLED NAVIER-STOKES PROBLEM----------------
//------------------------------------------------------------------------------
template<>
int Arlequin<2>::solveArlequinProblemMoving(int iterNumber, double tolerance,
                                            int problem_type, 
                                            int time_dependency){

    Mat               A,F;
    Vec               b, u, All;
    PetscErrorCode    ierr;
    PetscInt          Istart, Iend, Ii, Ione, iterations;
    KSP               ksp; 
    PC                pc;
    VecScatter        ctx;
    PetscScalar       val;
    //MatNullSpace      nullsp;
    PetscLogDouble bytes = 0;


    std::ofstream dragLift;
    dragLift.open("dragLift.dat", std::ofstream::out | std::ofstream::app);
    if (rank == 0) {
        dragLift << "Time   Pressure Drag   Pressure Lift " 
                 << "Friction Drag  Friction Lift Drag    Lift " 
                 << std::endl;
    };   

    if ((problem_type < 1) || (problem_type > 2)){
        std::cout << "WRONG PROBLEM TYPE." << std::endl;
        return 0;
    };

    if (time_dependency == 0) numTimeSteps = 1;
    
    iTimeStep = 0.;
    
    // //Construct the glue zone based on some defined criterion
    // setCouplingZone();

    // //Computes the Weight function for all the finite elements
    // setWeightFunction(16.);
   
    //Computes the Nodal correspondence between fine nodes and coarse elements
    // setNodalCorrespondenceFine();

    double integScheme = fineModel.integScheme;

    alpha_f = 1. / (1. + integScheme);
    alpha_m = 0.5 * (3. - integScheme) / (1. + integScheme);
    gamma = 0.5 + alpha_m - alpha_f;
    if (rank == 0) std::cout << "Time integ parameters: " << alpha_f << " " << alpha_m << " " << gamma << std::endl;

    // Computes the system size
    int sysSize = 3 * numNodesCoarse + 3 * numNodesFine + 2 * numNodesGlueZoneFine;
    
    numTimeSteps = 10000;

    for (iTimeStep = 0; iTimeStep < numTimeSteps; iTimeStep++){
        
        if (rank == 0) {std::cout << "------------------------- TIME STEP = "
                                  << iTimeStep << " -------------------------"
                                  << std::endl;}
        PetscMemoryGetCurrentUsage(&bytes);
        PetscPrintf(PETSC_COMM_WORLD,"Memory used %g M\n",bytes/(1024*1024));
        
        //Updates velocity and acceleration
        for (int i = 0; i < numNodesCoarse; i++){
            double accel[2], u[2], uprev[2];
            
            //Compute acceleration
            u[0] = nodesCoarse_[i] -> getVelocity(0);
            u[1] = nodesCoarse_[i] -> getVelocity(1);

            nodesCoarse_[i] -> setPreviousVelocity(u);
            
            accel[0] = nodesCoarse_[i] -> getAcceleration(0);
            accel[1] = nodesCoarse_[i] -> getAcceleration(1);
            
            nodesCoarse_[i] -> setPreviousAcceleration(accel);

            accel[0] *= (gamma - 1.) / gamma;
            accel[1] *= (gamma - 1.) / gamma;
            
            nodesCoarse_[i] -> setAcceleration(accel);            

        };

        // double f = .35;
        // double w = 2 * pi * f;

        for (int i = 0; i < numNodesFine; i++){
            double accel[2], u[2], uprev[2];
            
            //Compute acceleration
            u[0] = nodesFine_[i] -> getVelocity(0);
            u[1] = nodesFine_[i] -> getVelocity(1);

            nodesFine_[i] -> setPreviousVelocity(u);
            
            accel[0] = nodesFine_[i] -> getAcceleration(0);
            accel[1] = nodesFine_[i] -> getAcceleration(1);
            
            nodesFine_[i] -> setPreviousAcceleration(accel);

            accel[0] *= (gamma - 1.) / gamma;
            accel[1] *= (gamma - 1.) / gamma;
            
            nodesFine_[i] -> setAcceleration(accel);



            typename Nodes::VecLocD x, xn, xi;
            xi = nodesFine_[i] -> getInitialCoordinates();       
            x = nodesFine_[i] -> getCoordinates();       
    
            double a = -20 * pi / 180 + 10 * pi / 180 * cos(2.*pi*1.0*iTimeStep*dTime);// + 10 * pi / 180;

            // std::cout << " AAA " << a << std::endl;

            xn(0) = 0.5 + (xi(0)-0.5) * cos(a) - (xi(1)-0.0) * sin(a);
            xn(1) = 0.0 + (xi(0)-0.5) * sin(a) + (xi(1)-0.0) * cos(a);

            u[0] = (xn(0) - x(0)) / dTime;
            u[1] = (xn(1) - x(1)) / dTime;


            nodesFine_[i] -> setMeshVelocity(u);
      
            nodesFine_[i] -> setPreviousCoordinates(0,x(0));
            nodesFine_[i] -> setPreviousCoordinates(1,x(1));

            nodesFine_[i] -> setCoordinates(xn);
        };
        
        setSignaledDistance();
        setWeightFunction(1.);
        setNodalCorrespondenceFine();


        //STARTS NEWTON-RAPHSON
        for (int inewton = 0; inewton < iterNumber; inewton++){
            
            boost::posix_time::ptime t1 =
                boost::posix_time::microsec_clock::local_time();
            
            // Preallocates the matrix
            ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                                sysSize, sysSize, 300, NULL, 900, NULL, &A); 
            
            CHKERRQ(ierr);
            
            // Divides the matrix between the processes
            ierr = MatGetOwnershipRange(A, &Istart, &Iend);CHKERRQ(ierr);
            
            //Create PETSc vectors
            ierr = VecCreate(PETSC_COMM_WORLD, &b); CHKERRQ(ierr);
            ierr = VecSetSizes(b, PETSC_DECIDE, sysSize); CHKERRQ(ierr);
            ierr = VecSetFromOptions(b); CHKERRQ(ierr); 
            ierr = VecDuplicate(b, &u); CHKERRQ(ierr);
            ierr = VecDuplicate(b, &All); CHKERRQ(ierr);
                        
            for (int i=0; i<sysSize; i++){
                double val = 1.e-15;
                ierr = MatSetValues(A,1,&i,1,&i,&val,ADD_VALUES);
                
            }
                                

            //------------------------------------------------------------------
            //------------------BEGIN OF LINEAR SYSTEM ASSEMBLY-----------------
            //------------------------------------------------------------------

            //Coarse mesh
            
            for (int jel = 0; jel < numElemCoarse; jel++){   
                
                if (domDecompCoarse.first[jel] == rank) {
                    
                    //Compute Element matrix
                    std::pair<Elements::LocalMatrix, Elements::LocalVector> localMV;
                    localMV = elementsCoarse_[jel] -> getTransientNavierStokes();
        
                    typename Elements::Connectivity connec;
                    connec = elementsCoarse_[jel] -> getConnectivity();
                    
                    //Disperse local contributions into the global matrix
                    for (int i=0; i<6; i++){
                        for (int j=0; j<6; j++){
                            
                            if (fabs(localMV.first(2*i  ,2*j  )) >= 1.e-15){
                                int dof_i = 2 * connec(i);
                                int dof_j = 2 * connec(j);
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&localMV.first(2*i  ,2*j  ),ADD_VALUES);
                            };
                            if (fabs(localMV.first(2*i+1,2*j  )) >= 1.e-15){
                                int dof_i = 2 * connec(i) + 1;
                                int dof_j = 2 * connec(j);
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&localMV.first(2*i+1,2*j  ),ADD_VALUES);
                            };

                            if (fabs(localMV.first(2*i  ,2*j+1)) >= 1.e-15){
                                int dof_i = 2 * connec(i);
                                int dof_j = 2 * connec(j) + 1;
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&localMV.first(2*i  ,2*j+1),ADD_VALUES);
                            };

                            if (fabs(localMV.first(2*i+1,2*j+1)) >= 1.e-15){
                                int dof_i = 2 * connec(i) + 1;
                                int dof_j = 2 * connec(j) + 1;
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&localMV.first(2*i+1,2*j+1),ADD_VALUES);
                            };

                            //Matrix Q and Qt
                            if (fabs(localMV.first(2*i  ,12+j)) >= 1.e-15){
                                int dof_i = 2 * connec(i);
                                int dof_j = 2 * numNodesCoarse + connec(j);
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&localMV.first(2*i  ,12+j),ADD_VALUES);
                            };
                            
                            if (fabs(localMV.first(12+j,2*i  )) >= 1.e-15){
                                int dof_i = 2 * connec(i);
                                int dof_j = 2 * numNodesCoarse + connec(j);
                                ierr = MatSetValues(A,1,&dof_j,1,&dof_i,&localMV.first(12+j,2*i  ),ADD_VALUES);
                            };

                            if (fabs(localMV.first(2*i+1,12+j)) >= 1.e-15){
                                int dof_i = 2 * connec(i) + 1;
                                int dof_j = 2 * numNodesCoarse + connec(j);
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&localMV.first(2*i+1,12+j),ADD_VALUES);
                            };

                            if (fabs(localMV.first(12+j,2*i+1)) >= 1.e-15){
                                int dof_i = 2 * connec(i) + 1;
                                int dof_j = 2 * numNodesCoarse + connec(j);
                                ierr = MatSetValues(A,1,&dof_j,1,&dof_i,&localMV.first(12+j,2*i+1),ADD_VALUES);
                            };

                            if (fabs(localMV.first(12+i,12+j)) >= 1.e-15){
                                int dof_i = 2 * numNodesCoarse + connec(i);
                                int dof_j = 2 * numNodesCoarse + connec(j);
                                ierr = MatSetValues(A,1,&dof_j,1,&dof_i,&localMV.first(12+i,12+j),ADD_VALUES);
                            };
                        };                  
                        //Rhs vector
                        if (fabs(localMV.second(2*i  )) >= 1.e-15){
                            int dof_i = 2 * connec(i);
                            ierr = VecSetValues(b,1,&dof_i,&localMV.second(2*i  ),ADD_VALUES);
                        };

                        if (fabs(localMV.second(2*i+1)) >= 1.e-15){
                            int dof_i = 2 * connec(i) + 1;
                            ierr = VecSetValues(b,1,&dof_i,&localMV.second(2*i+1),ADD_VALUES);
                        };

                        if (fabs(localMV.second(12+i)) >= 1.e-15){
                            int dof_i = 2 * numNodesCoarse + connec(i);
                            ierr = VecSetValues(b,1,&dof_i,&localMV.second(12+i),ADD_VALUES);
                        };
                    };
                };
            };

            //Fine mesh
            for (int jel = 0; jel < numElemFine; jel++){   
                
                if (domDecompFine.first[jel] == rank) {
                    
                    //Compute Element matrix
                    std::pair<Elements::LocalMatrix, Elements::LocalVector> localMV;
                    localMV = elementsFine_[jel] -> getTransientNavierStokes();
        
                    typename Elements::Connectivity connec;
                    connec = elementsFine_[jel] -> getConnectivity();          

                    //Disperse local contributions into the global matrix
                    for (int i=0; i<6; i++){
                        for (int j=0; j<6; j++){
                            if (fabs(localMV.first(2*i  ,2*j  )) >= 1.e-15){
                                int dof_i = 3*numNodesCoarse + 2 * connec(i);
                                int dof_j = 3*numNodesCoarse + 2 * connec(j);
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&localMV.first(2*i  ,2*j  ),ADD_VALUES);
                            };
                            
                            if (fabs(localMV.first(2*i+1,2*j  )) >= 1.e-15){
                                int dof_i = 3*numNodesCoarse + 2 * connec(i) +1;
                                int dof_j = 3*numNodesCoarse + 2 * connec(j);
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&localMV.first(2*i+1,2*j  ),ADD_VALUES);
                            };

                            if (fabs(localMV.first(2*i  ,2*j+1)) >= 1.e-15){
                                int dof_i = 3*numNodesCoarse + 2 * connec(i);
                                int dof_j = 3*numNodesCoarse + 2 * connec(j) +1;
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&localMV.first(2*i  ,2*j+1),ADD_VALUES);
                            };

                            if (fabs(localMV.first(2*i+1,2*j+1)) >= 1.e-15){
                                int dof_i = 3*numNodesCoarse + 2 * connec(i) +1;
                                int dof_j = 3*numNodesCoarse + 2 * connec(j) +1;
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&localMV.first(2*i+1,2*j+1),ADD_VALUES);
                            };

                            //Matrix Q and Qt
                            if (fabs(localMV.first(2*i  ,12+j)) >= 1.e-15){
                                int dof_i = 3 * numNodesCoarse + 2 * connec(i);
                                int dof_j = 3 * numNodesCoarse + 2 * numNodesFine + connec(j);
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&localMV.first(2*i  ,12+j),ADD_VALUES);
                            };

                            if (fabs(localMV.first(12+j,2*i  )) >= 1.e-15){
                                int dof_i = 3 * numNodesCoarse + 2 * connec(i);
                                int dof_j = 3 * numNodesCoarse + 2 * numNodesFine + connec(j);
                                ierr = MatSetValues(A,1,&dof_j,1,&dof_i,&localMV.first(12+j,2*i  ),ADD_VALUES);
                            };

                            if (fabs(localMV.first(2*i+1,12+j)) >= 1.e-15){
                                int dof_i = 3 * numNodesCoarse + 2 * connec(i) + 1;
                                int dof_j = 3 * numNodesCoarse + 2 * numNodesFine + connec(j);
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&localMV.first(2*i+1,12+j),ADD_VALUES);
                            };
                            
                            if (fabs(localMV.first(12+j,2*i+1)) >= 1.e-15){
                                int dof_i = 3 * numNodesCoarse + 2 * connec(i) + 1;
                                int dof_j = 3 * numNodesCoarse + 2 * numNodesFine + connec(j);
                                ierr = MatSetValues(A,1,&dof_j,1,&dof_i,&localMV.first(12+j,2*i+1),ADD_VALUES);
                            };
                            
                            if (fabs(localMV.first(12+i,12+j)) >= 1.e-15){
                                int dof_i = 3 * numNodesCoarse + 2 * numNodesFine + connec(i);
                                int dof_j = 3 * numNodesCoarse + 2 * numNodesFine + connec(j);
                                ierr = MatSetValues(A,1,&dof_j,1,&dof_i,&localMV.first(12+i,12+j),ADD_VALUES);
                            };
                        };     
                        ///Rhs vector
                        if (fabs(localMV.second(2*i  )) >= 1.e-15){
                            int dof_i = 3 * numNodesCoarse + 2 * connec(i);
                            ierr = VecSetValues(b,1,&dof_i,&localMV.second(2*i  ),ADD_VALUES);
                        };

                        if (fabs(localMV.second(2*i+1)) >= 1.e-15){
                            int dof_i = 3 * numNodesCoarse + 2 * connec(i) + 1;
                            ierr = VecSetValues(b,1,&dof_i,&localMV.second(2*i+1),ADD_VALUES);
                        };

                        if (fabs(localMV.second(12+i)) >= 1.e-15){
                            int dof_i = 3 * numNodesCoarse + 2 * numNodesFine + connec(i);
                            ierr = VecSetValues(b,1,&dof_i,&localMV.second(12+i),ADD_VALUES);
                        };
                    }; 
                };                
            };
              

            //Lagrange Multipliers
            for (int l=0; l< numElemGlueZoneFine; l++){
                
                int jel = elementsGlueZoneFine_[l];
                
                if (domDecompFine.first[jel] == rank) {
                    
                    typename Elements::LocalMatrix Ajac, AjacAnt, AStab, ArlequinM, ArlequinM2;
                    typename Elements::LocalVector Rhs, RhsStab, RhsArlequin,rhsLagMult;
                    typename Elements::Connectivity connec, connecC, connecL;
                    
                    connec = elementsFine_[jel] -> getConnectivity();
                    connecL = glueZoneFine_[l] -> getConnectivity();

                    // FINE MESH
                    
                    //Computes element matrix
                    elementsFine_[jel] -> getLagrangeMultipliersSameMesh();
                    
                    //Gets element matrice and rhs vectors
                    // -L1
                    Ajac = - elementsFine_[jel] -> getLagrMultMatrix();
                    // +L1 * Lambda
                    // rhsLagMult = - elementsFine_[jel] -> getRhsLagrangeMultipliers(-Ajac);
                    rhsLagMult = - elementsFine_[jel] -> getLagrMultVector();

                    //PSPG and SUPG stabilizations
                    AStab = elementsFine_[jel] -> getJacNRMatrix();
                    RhsStab = elementsFine_[jel] -> getRhsVector();

                    //Arlequin Stabilization
                    ArlequinM = elementsFine_[jel] -> getArlequinStabilizationMatrix();
                    RhsArlequin = elementsFine_[jel] -> getArlequinStabilizationVector();

                    ArlequinM2 = elementsFine_[jel] -> getArlequinStabilizationMatrix2();

                    // -L1t
                    Ajac = trans(Ajac);

                    // L1t * u1
                    Rhs = -elementsFine_[jel] -> getRhsVelocities(-Ajac);

                    std::pair<Elements::LocalVector, Elements::LocalMatrix>  lagMult;
   
                    lagMult.first.clear();
                    lagMult.second.clear();
                         
                    lagMult = elementsFine_[jel] -> getBoundaryConditionsVelocity(Rhs,Ajac);

                    Ajac = lagMult.second;
                    Rhs = lagMult.first;


                    lagMult.first.clear();
                    lagMult.second.clear();
                         
                    lagMult = elementsFine_[jel] -> getBoundaryConditionsLagMult(rhsLagMult,Ajac);

                    Ajac = lagMult.second;
                    rhsLagMult = lagMult.first;


                    double integ = alpha_f * gamma * dTime;
                    //Disperse local contributions into the global matrix
                    for (int i=0; i<6; i++){
                        for (int j=0; j<6; j++){
                            
                            //COUPLING OPERATOR
                            if (fabs(Ajac(2*i  ,2*j  )) >= 1.e-15){
                                int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                                int d_j = 3*numNodesCoarse + 2*connec(j);
                                double value = Ajac(2*i  ,2*j  )*integ;
                                ierr = MatSetValues(A,1,&d_i,1,&d_j,&value,ADD_VALUES);
                                ierr = MatSetValues(A,1,&d_j,1,&d_i,&Ajac(2*i  ,2*j  ),ADD_VALUES);
                            };
                            if (fabs(Ajac(2*i+1,2*j  )) >= 1.e-15){
                                int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i) + 1;
                                int d_j = 3*numNodesCoarse + 2*connec(j);
                                double value = Ajac(2*i+1,2*j  ) * integ;
                                ierr = MatSetValues(A,1,&d_i,1,&d_j,&value,ADD_VALUES);
                                ierr = MatSetValues(A,1,&d_j,1,&d_i,&Ajac(2*i+1,2*j  ),ADD_VALUES);
                            };
                            if (fabs(Ajac(2*i+1,2*j+1)) >= 1.e-15){
                                int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i) + 1;
                                int d_j = 3*numNodesCoarse + 2*connec(j) + 1;
                                double value = Ajac(2*i+1,2*j+1)*integ;
                                ierr = MatSetValues(A,1,&d_i,1,&d_j,&value,ADD_VALUES);
                                ierr = MatSetValues(A,1,&d_j,1,&d_i,&Ajac(2*i+1,2*j+1),ADD_VALUES);
                            };
                            if (fabs(Ajac(2*i  ,2*j+1)) >= 1.e-15){
                                int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                                int d_j = 3 * numNodesCoarse + 2*connec(j) + 1;
                                double value = Ajac(2*i  ,2*j+1)*integ;
                                ierr = MatSetValues(A,1,&d_i,1,&d_j,&value,ADD_VALUES);
                                ierr = MatSetValues(A,1,&d_j,1,&d_i,&Ajac(2*i  ,2*j+1),ADD_VALUES);
                            };

                            //SUPG STABILIZATION
                            if (fabs(AStab(2*i  ,2*j  )) >= 1.e-15){
                                int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                                int d_j = 3*numNodesCoarse + 2*connec(j);
                                ierr = MatSetValues(A,1,&d_j,1,&d_i,&AStab(2*i  ,2*j  ),ADD_VALUES);
                            };
                            if (fabs(AStab(2*i+1,2*j+1)) >= 1.e-15){
                                int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i) + 1;
                                int d_j = 3*numNodesCoarse + 2*connec(j) + 1;
                                ierr = MatSetValues(A,1,&d_j,1,&d_i,&AStab(2*i+1,2*j+1),ADD_VALUES);
                            };

                            //PSPG STABILIZATION
                            if (fabs(AStab(2*i  ,12+j)) >= 1.e-15){
                                int dof_i = 3*numNodesCoarse + 2*numNodesFine + connec(i);
                                int dof_j = 3*numNodesCoarse + 3*numNodesFine + connecL(j);
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&AStab(2*i  ,12+j),ADD_VALUES);
                            };
                            if (fabs(AStab(2*i+1,12+j)) >= 1.e-15){
                                int dof_i = 3*numNodesCoarse + 2*numNodesFine + connec(i);
                                int dof_j = 3*numNodesCoarse + 3*numNodesFine + connecL(j)+1;
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&AStab(2*i+1,12+j),ADD_VALUES);
                            };

                            //ARLEQUIN STABILIZATION
                            // if (fabs(ArlequinM(2*i  ,2*j+1)) >= 1.e-15){
                            //     int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                            //     int d_j = 3*numNodesCoarse + 2*connec(j) + 1;
                            //     ierr = MatSetValues(A,1,&d_i,1,&d_j,
                            //                         &ArlequinM(2*i  ,2*j+1),ADD_VALUES);
                            // };
                            if (fabs(ArlequinM2(2*i  ,2*j  )) >= 1.e-15){
                                int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                                int d_j = 3*numNodesCoarse + 2*connec(j);
                                ierr = MatSetValues(A,1,&d_i,1,&d_j,&ArlequinM2(2*i  ,2*j  ),ADD_VALUES);
                            };
                            // if (fabs(ArlequinM(2*i+1,2*j  )) >= 1.e-15){
                            //     int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i) + 1;
                            //     int d_j = 3*numNodesCoarse + 2*connec(j);
                            //     ierr = MatSetValues(A,1,&d_i,1,&d_j,
                            //                         &ArlequinM(2*i+1,2*j  ),ADD_VALUES);
                            // };
                            if (fabs(ArlequinM2(2*i+1,2*j+1)) >= 1.e-15){
                                int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i) + 1;
                                int d_j = 3*numNodesCoarse + 2*connec(j) + 1;
                                ierr = MatSetValues(A,1,&d_i,1,&d_j,&ArlequinM2(2*i+1,2*j+1),ADD_VALUES);
                            };

                            if (fabs(ArlequinM2(12+i,2*j  )) >= 1.e-15){
                                int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                                int dof_j = 3*numNodesCoarse + 2*numNodesFine + connec(j);
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinM2(12+i,2*j  ),ADD_VALUES);
                            };
                            if (fabs(ArlequinM2(12+i,2*j+1)) >= 1.e-15){
                                int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i) + 1;
                                int dof_j = 3*numNodesCoarse + 2*numNodesFine + connec(j);
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinM2(12+i,2*j+1),ADD_VALUES);
                            };
                            if (fabs(ArlequinM(2*i  ,2*j  )) >= 1.e-15){
                                int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                                int dof_j = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(j);
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinM(2*i  ,2*j  ),ADD_VALUES);
                                dof_i++; dof_j++;
                                ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinM(2*i+1,2*j+1),ADD_VALUES);
                            };                            
                        };

                        //RHS VECTOR
                        //COUPLING OPERATOR
                        int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                        ierr = VecSetValues(b,1,&dof_i,&Rhs(2*i  ),ADD_VALUES);

                        dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i) + 1;
                        ierr = VecSetValues(b,1,&dof_i,&Rhs(2*i+1),ADD_VALUES);

                        dof_i = 3*numNodesCoarse + 2*connec(i);
                        ierr = VecSetValues(b,1,&dof_i,&rhsLagMult(2*i  ),ADD_VALUES);

                        dof_i = 3*numNodesCoarse + 2*connec(i) + 1;
                        ierr = VecSetValues(b,1,&dof_i,&rhsLagMult(2*i+1),ADD_VALUES);

                        //SUPG STABILIZATION
                        dof_i = 3*numNodesCoarse + 2*connec(i);                       
                        ierr = VecSetValues(b,1,&dof_i,&RhsStab(2*i  ),ADD_VALUES);

                        dof_i = 3*numNodesCoarse + 2*connec(i) + 1;
                        ierr = VecSetValues(b,1,&dof_i,&RhsStab(2*i+1),ADD_VALUES);

                        //PSPG STABILIZATION
                        dof_i = 3*numNodesCoarse + 2*numNodesFine + connec(i);
                        ierr = VecSetValues(b,1,&dof_i,&RhsStab(12+i),ADD_VALUES);

                        // ///ARLEQUIN STABILIZATION
                        dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                        ierr = VecSetValues(b,1,&dof_i,&RhsArlequin(2*i  ),ADD_VALUES);
                        dof_i++;
                        ierr = VecSetValues(b,1,&dof_i,&RhsArlequin(2*i+1),ADD_VALUES);
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
                        double pspg = 0.;//elementsCoarse_[iElemCoarse] -> getPSPG();
                        ublas::bounded_vector<double, 6> press_, velX_, velY_;


                        for (int k = 0; k < 6; k++){
                            press_(k) = nodesCoarse_[connecC(k)] -> getPressure();
                            velX_(k) = nodesCoarse_[connecC(k)] -> getVelocity(0);
                            velY_(k) = nodesCoarse_[connecC(k)] -> getVelocity(1);
                        }
                        
                        elementsFine_[jel] -> 
                            getLagrangeMultipliersDifferentMesh(iElemCoarse,pspg,press_,velX_,velY_);
                        
                        //Computes element matrix
                        // L0
                        Ajac = elementsFine_[jel] -> getLagrMultMatrix();
                        // L0t
                        Ajac = trans(Ajac);
                        AStab = elementsFine_[jel] -> getJacNRMatrix();
                        AStab = trans(AStab);
                        RhsStab = elementsFine_[jel] -> getRhsVector();
                        ArlequinM = elementsFine_[jel] -> getArlequinStabilizationMatrix();
                        ArlequinM = trans(ArlequinM);
                        RhsArlequin = elementsFine_[jel] -> getArlequinStabilizationVector();
                        
                        ArlequinM2 = elementsFine_[jel] -> getArlequinStabilizationMatrix2();
                        ArlequinM2 = trans(ArlequinM2);

                        // -L0 * lambda
                        rhsLagMult =  elementsFine_[jel] -> getRhsLagrangeMultipliers(trans(Ajac));

                        connecC = elementsCoarse_[iElemCoarse] -> 
                            getConnectivity();

                        // AStab.clear();
                        // // RhsStab.clear();
                        // ArlequinM.clear();
                        // RhsArlequin.clear();


                        
                        // U_.clear();
                        // for (int i = 0; i < 6; i++){
                        //     U_(2*i  ) = nodesFine_[connec(i)] -> getLagrangeMultiplier(0);
                        //     U_(2*i+1) = nodesFine_[connec(i)] -> getLagrangeMultiplier(1);
                        // };
                        
                        // noalias(rhsLagMult) = -prod(trans(Ajac),U_);
                        // noalias(lagStab) = -prod(trans(AjacAnt),U_);
                        


 
                    Rhs = elementsCoarse_[iElemCoarse] -> getRhsVelocities(Ajac);

                    std::pair<Elements::LocalVector, Elements::LocalMatrix>  lagMult;
                    


                    // lagMult.first.clear();
                    // lagMult.second.clear();
                         
                    // lagMult = elementsCoarse_[iElemCoarse] -> getBoundaryConditionsVelocity(rhsLagMult,(Ajac));

                    // Ajac = lagMult.second;
                    // rhsLagMult = lagMult.first;


                    // lagMult.first.clear();
                    // lagMult.second.clear();
                         
                    // lagMult = elementsFine_[jel] -> getBoundaryConditionsLagMult(Rhs,(Ajac));

                    // Ajac = (lagMult.second);
                    // Rhs = lagMult.first;


                    // if (jel == 94){
                    //     std::cout << "MATRIZ " << std::endl;
                    //     for (int i=0; i<18; i++){
                    //         for (int j=0; j<18; j++){
                    //             std::cout << Ajac(i,j) << " ";
                    //         }
                    //         std::cout << std::endl;
                    //     }
                    //     std::cout << "VETOR " << std::endl;
                    //     for (int i=0; i<18; i++){
                    //         std::cout << Rhs(i) << " " << rhsLagMult(i) << std::endl;
                    //     }
                    // }

                    rhsLagMult *= 1000.e0;
                    //Ajac *= 1000.;

                    for (int i = 0; i < 12; i++){
                        RhsStab(i) *= 1000.e0;
                        for (int j = 0; j < 12; j++){
                            AStab(i,j)*=1000.e0;
                            Ajac(i,j)*=1000.e0;
                            ArlequinM2(i,j)*=1000.e0;
                        }
                    }
                        // std::pair<Elements::LocalVector, 
                        //           Elements::LocalMatrix>  lagMult;
                        
                        // lagMult.first.clear();
                        // lagMult.second.clear();
                         
                        // lagMult = elementsCoarse_[iElemCoarse] -> 
                        //     getRhsVectorAndBoundaryConditions(Ajac);

                        // noalias(rhsLagMult) = -prod((lagMult.second),U_);
                    // for (int i=0; i<6; i++){
                    //     if(lagMult.second(12+i,12+i) > 0){
                    //         rhsLagMult(2*i  ) = 0.;
                    //         rhsLagMult(2*i+1) = 0.;
                    //         double one = 1.;
                    //         int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                    //         ierr = MatSetValues(A,1,&d_i,1,&d_i,&one,INSERT_VALUES);
                    //         d_i++;
                    //         ierr = MatSetValues(A,1,&d_i,1,&d_i,&one,INSERT_VALUES);
                    //     }
                    // }
                        double integ = alpha_f * gamma * dTime;
                        //Disperse local contribution into the global matrix
                        for (int i=0; i<6; i++){
                            for (int j=0; j<6; j++){
                                //COUPLING OPERATOR
                                if (fabs(Ajac(2*i  ,2*j  )) >= 1.e-15){
                                    int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                                    int dof_j = 2*connecC(j);
                                    double value = Ajac(2*i  ,2*j  ) * integ;
                                    ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&value,ADD_VALUES);
                                    ierr = MatSetValues(A,1,&dof_j,1,&dof_i,&Ajac(2*i  ,2*j  ),ADD_VALUES);
                                }
                                if (fabs(Ajac(2*i+1,2*j  )) >= 1.e-15){
                                    int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i) + 1;
                                    int dof_j = 2*connecC(j);
                                    double value = Ajac(2*i+1,2*j  )*integ;
                                    ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&value,ADD_VALUES);
                                    ierr = MatSetValues(A,1,&dof_j,1,&dof_i,&Ajac(2*i+1,2*j  ),ADD_VALUES);
                                };
                                if (fabs(Ajac(2*i+1,2*j+1)) >= 1.e-15){
                                    int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i) + 1;
                                    int dof_j = 2*connecC(j) + 1;
                                    double value = Ajac(2*i+1,2*j+1)*integ;
                                    ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&value,ADD_VALUES);
                                    ierr = MatSetValues(A,1,&dof_j,1,&dof_i,&Ajac(2*i+1,2*j+1),ADD_VALUES);
                                };
                                if (fabs(Ajac(2*i  ,2*j+1)) >= 1.e-15){
                                    int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                                    int dof_j = 2*connecC(j) + 1;
                                    double value = Ajac(2*i  ,2*j+1)*integ;
                                    ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&value,ADD_VALUES);
                                    ierr = MatSetValues(A,1,&dof_j,1,&dof_i,&Ajac(2*i  ,2*j+1),ADD_VALUES);
                                };
                                if (fabs(Ajac(12+i,12+j)) >= 1.e-15){
                                    int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                                    int dof_j = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(j);
                                    // ierr = MatSetValues(A,1,&dof_i,1,&dof_j,
                                    //                 &lagMult.second(12+i,12+j),ADD_VALUES);
                                    dof_i++; dof_j++;
                                    ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&Ajac(12+i,12+j),ADD_VALUES);
                                };
                         

                                //SUPG STABILIZATION
                                if (fabs(AStab(2*i  ,2*j  )) >= 1.e-15){
                                    int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                                    int dof_j = 2*connecC(j);
                                    ierr = MatSetValues(A,1,&dof_j,1,&dof_i,&AStab(2*i  ,2*j  ),ADD_VALUES);
                                };
                                if (fabs(AStab(2*i+1,2*j+1)) >= 1.e-15){
                                    int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i) + 1;
                                    int dof_j = 2*connecC(j) + 1;
                                    ierr = MatSetValues(A,1,&dof_j,1,&dof_i,&AStab(2*i+1,2*j+1),ADD_VALUES);
                                };

                                //PSPG STABILIZATION
                                if (fabs(AStab(2*i  ,12+j)) >= 1.e-15){
                                    int dof_i = 2*numNodesCoarse + connecC(i);
                                    int dof_j = 3*numNodesCoarse + 2*numNodesFine + connecL(j);
                                    ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&AStab(2*i  ,12+j),ADD_VALUES);
                                };
                                if (fabs(AStab(2*i+1,12+j)) >= 1.e-15){
                                    int dof_i = 2*numNodesCoarse + connecC(i);
                                    int dof_j = 3*numNodesCoarse + 2*numNodesFine + connecL(j)+1;
                                    ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&AStab(2*i+1,12+j),ADD_VALUES);
                                };

                                ///ARLEQUIN STABILIZATION
                                // if (fabs(ArlequinM(2*i  ,2*j+1)) >= 1.e-15){
                                //     int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                                //     int d_j = 2*connecC(j) + 1;
                                //     ierr = MatSetValues(A,1,&d_i,1,&d_j,
                                //                         &ArlequinM(2*i  ,2*j+1),ADD_VALUES);
                                // };
                                if (fabs(ArlequinM2(2*i  ,2*j  )) >= 1.e-15){
                                    int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                                    int d_j = 2*connecC(j);
                                    ierr = MatSetValues(A,1,&d_i,1,&d_j,&ArlequinM2(2*i  ,2*j  ),ADD_VALUES);
                                };
                                // if (fabs(ArlequinM(2*i+1,2*j  )) >= 1.e-15){
                                //     int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i) + 1;
                                //     int d_j = 2*connecC(j);
                                //     ierr = MatSetValues(A,1,&d_i,1,&d_j,
                                //                         &ArlequinM(2*i+1,2*j  ),ADD_VALUES);
                                // };
                                if (fabs(ArlequinM2(2*i+1,2*j+1)) >= 1.e-15){
                                    int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i) + 1;
                                    int d_j = 2*connecC(j) + 1;
                                    ierr = MatSetValues(A,1,&d_i,1,&d_j,&ArlequinM2(2*i+1,2*j+1),ADD_VALUES);
                                };

                                if (fabs(ArlequinM2(12+i,2*j  )) >= 1.e-15){
                                    int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                                    int dof_j = 2*numNodesCoarse + connecC(j);
                                    ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinM2(12+i,2*j  ),ADD_VALUES);
                                };
                                if (fabs(ArlequinM2(12+i,2*j+1)) >= 1.e-15){
                                    int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i) + 1;
                                    int dof_j = 2*numNodesCoarse + connecC(j);
                                    ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinM2(12+i,2*j+1),ADD_VALUES);
                                };

                                if (fabs(ArlequinM(2*i  ,2*j  )) >= 1.e-15){
                                    int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                                    int dof_j = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(j);
                                    ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinM(2*i  ,2*j  ),ADD_VALUES);
                                    dof_i++; dof_j++;
                                    ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinM(2*i+1,2*j+1),ADD_VALUES);
                                };
                            };
                            //RHS VECTOR
                            //COUPLING OPERATOR
                            int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                            ierr = VecSetValues(b,1,&d_i,&Rhs(2*i  ),ADD_VALUES);

                            d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i) + 1;
                            ierr = VecSetValues(b,1,&d_i,&Rhs(2*i+1),ADD_VALUES);

                            int dof_i = 2*connecC(i);
                            ierr = VecSetValues(b,1,&dof_i,&rhsLagMult(2*i  ),ADD_VALUES);

                            dof_i = 2*connecC(i) + 1;
                            ierr = VecSetValues(b,1,&dof_i,&rhsLagMult(2*i+1),ADD_VALUES);

                            //SUPG STABILIZATION
                            dof_i = 2*connecC(i);
                            ierr = VecSetValues(b,1,&dof_i,&RhsStab(2*i  ),ADD_VALUES);

                            dof_i = 2*connecC(i) + 1;
                            ierr = VecSetValues(b,1,&dof_i,&RhsStab(2*i+1),ADD_VALUES);

                            //PSPG STABILIZATION
                            dof_i = 2*numNodesCoarse + connecC(i);
                            ierr = VecSetValues(b,1,&dof_i,&RhsStab(12+i),ADD_VALUES);

                            //ARLEQUIN STABILIZATION
                            dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL(i);
                            ierr = VecSetValues(b,1,&dof_i,&RhsArlequin(2*i  ),ADD_VALUES);
                            dof_i++;
                            ierr = VecSetValues(b,1,&dof_i,&RhsArlequin(2*i+1),ADD_VALUES);
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

        ierr = PCFactorSetMatSolverType(pc,MATSOLVERMUMPS);
        PCFactorSetUpMatSolverType(pc);
        PCFactorGetMatrix(pc,&F);

        PetscInt ival,icntl;
        icntl = 14; ival = 60;
        MatMumpsSetIcntl(F,icntl,ival);
        
#endif
        
        
        ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
        ierr = KSPSetUp(ksp);
        
        
#if defined(PETSC_HAVE_MUMPS)
        PetscInt  info1,info2,icntl14;

        MatMumpsGetInfo(F,1,&info1);
        MatMumpsGetInfo(F,2,&info2);
        MatMumpsGetIcntl(F,14,&icntl14);
        if((rank == 0) && (info1 != 0)) std::cout << " INFO(1) = " << info1
                                                      << " " << info2 << " " 
                                                      << icntl14 << std::endl;
#endif
            
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
            double normU = 0.;
            double normP = 0.;
            double normL = 0.;
            double normT = 0.;
            double p_;
            Ione = 1;


            for (int i = 0; i < numNodesCoarse; ++i){
                Ii = 2 * i;
                ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                if (nodesCoarse_[i] -> getDistFunction() > -1.2) val *= 1000.e0;
                u_[0] = val;
                normU += val*val;
                nodesCoarse_[i] -> incrementAcceleration(0,u_[0]);
                nodesCoarse_[i] -> incrementVelocity(0,u_[0]*gamma*dTime);
                
                Ii = 2 * i + 1;
                ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                if (nodesCoarse_[i] -> getDistFunction() > -1.2) val *= 1000.e0;
                u_[1] = val;
                normU += val*val;
                nodesCoarse_[i] -> incrementAcceleration(1,u_[1]);
                nodesCoarse_[i] -> incrementVelocity(1,u_[1]*gamma*dTime);
            };
            
            for (int i = 0; i<numNodesCoarse; i++){
                Ii = 2 * numNodesCoarse + i;
                ierr = VecGetValues(All,Ione,&Ii,&val);CHKERRQ(ierr);
                p_ = val;
                normP += val*val;
                nodesCoarse_[i] -> incrementPressure(p_);
            };
            
            for (int i = 0; i < numNodesFine; ++i){
                Ii = 3 * numNodesCoarse + 2 * i;
                ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                u_[0] = val;
                normU += val*val;
                nodesFine_[i] -> incrementAcceleration(0,u_[0]);
                nodesFine_[i] -> incrementVelocity(0,u_[0]*gamma*dTime);
                
                Ii = 3 * numNodesCoarse + 2 * i + 1;
                ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                u_[1] = val;
                normU += val*val;
                nodesFine_[i] -> incrementAcceleration(1,u_[1]);
                nodesFine_[i] -> incrementVelocity(1,u_[1]*gamma*dTime);
            };
            
            for (int i = 0; i<numNodesFine; i++){
                Ii = 3 * numNodesCoarse + 2 * numNodesFine + i;
                ierr = VecGetValues(All,Ione,&Ii,&val);CHKERRQ(ierr);
                p_ = val;
                normP += val*val;
                nodesFine_[i] -> incrementPressure(p_);
            };
            
            for (int i = 0; i < numNodesGlueZoneFine; ++i){
                Ii = 3 * numNodesCoarse + 3 * numNodesFine + 2 * i;
                ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                u_[0] = val;
                nodesFine_[nodesGlueZoneFine_[i]] -> 
                    incrementLagrangeMultiplier(0,u_[0]);
                normL += val*val;

                Ii = 3 * numNodesCoarse + 3 * numNodesFine + 2 * i + 1;
                ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                u_[1] = val;
                nodesFine_[nodesGlueZoneFine_[i]] -> 
                    incrementLagrangeMultiplier(1,u_[1]);
                normL += val*val;
                // std::cout << "LAG M " << u_[0] << " " << u_[1] << std::endl;
            };


            //Computes the solution vector norm
            ierr = VecNorm(u,NORM_2,&val);CHKERRQ(ierr);
            
            boost::posix_time::ptime t2 =
                boost::posix_time::microsec_clock::local_time();
            
            if(rank == 0){
                boost::posix_time::time_duration diff = t2 - t1;
                
                std::cout<<"Iteration = " << inewton << " (" << iterations <<  
                    ")  Du Norm = " << std::scientific << sqrt(normU) 
                         << " " << sqrt(normP) 
                         << " " << sqrt(normL) 
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
            nodesFine_[i] -> setPressureArlequin(nodesFine_[i] ->getPressure());
        };
        
        for (int i = 0; i<numNodesGlueZoneFine; i++){
            typename Nodes::VecLocD xsi;
            typename Quadrature::NodalValuesQuad u_coarse, v_coarse, p_coarse;
            typename Elements::Connectivity connecCoarse;
            
            double u = 0.;
            double v = 0.;
            double p = 0.;
            
            int elCoarse = nodesFine_[nodesGlueZoneFine_[i]] -> 
                getNodalElemCorrespondence();
            xsi = nodesFine_[nodesGlueZoneFine_[i]] -> 
                getNodalXsiCorrespondence();
            
            connecCoarse = elementsCoarse_[elCoarse] -> getConnectivity();
            
            for (int j=0; j<6; j++){
                u_coarse(j) = nodesCoarse_[connecCoarse(j)] -> getVelocity(0);
                v_coarse(j) = nodesCoarse_[connecCoarse(j)] -> getVelocity(1);
                p_coarse(j) = nodesCoarse_[connecCoarse(j)] -> getPressure();
            };
            
            shapeQuad.evaluate(xsi,phi_);
            
            for (int j=0; j<6; j++){
                u += u_coarse(j) * phi_(j);
                v += v_coarse(j) * phi_(j);
                p += p_coarse(j) * phi_(j);
            };
            
            double wFunc = nodesFine_[nodesGlueZoneFine_[i]] -> 
                getWeightFunction();
            
            double u_int = nodesFine_[nodesGlueZoneFine_[i]] -> 
                getVelocity(0) * wFunc + u * (1. - wFunc);
            double v_int = nodesFine_[nodesGlueZoneFine_[i]] -> 
                getVelocity(1) * wFunc + v * (1. - wFunc);
            double p_int = nodesFine_[nodesGlueZoneFine_[i]] -> 
                getPressure() * wFunc + p * (1. - wFunc);
            
            nodesFine_[nodesGlueZoneFine_[i]] -> setVelocityArlequin(0,u_int);
            nodesFine_[nodesGlueZoneFine_[i]] -> setVelocityArlequin(1,v_int);
            nodesFine_[nodesGlueZoneFine_[i]] -> setPressureArlequin(p_int);
            
        };
        
        for (int i=0; i<numNodesCoarse; i++){
            nodesCoarse_[i] -> 
                setVelocityArlequin(0,nodesCoarse_[i] -> getVelocity(0));
            nodesCoarse_[i] -> 
                setVelocityArlequin(1,nodesCoarse_[i] -> getVelocity(1));
            nodesCoarse_[i] -> 
                setPressureArlequin(nodesCoarse_[i] -> getPressure());
        };
        
        // Compute and print drag and lift coefficients
        if (fineModel.getComputeDragAndLift()){
            dragAndLiftCoefficients(dragLift);
        };

        if (rank == 0) {
            if (coarseModel.printVorticity){
                for (int i = 0; i < numNodesCoarse; i++){
                    nodesCoarse_[i] -> clearVorticity();
                };
                for (int i = 0; i < numNodesFine; i++){
                    nodesFine_[i] -> clearVorticity();
                };
                for (int jel = 0; jel < numElemCoarse; jel++){
                    elementsCoarse_[jel] -> computeVorticity();
                };
                for (int jel = 0; jel < numElemFine; jel++){
                    elementsFine_[jel] -> computeVorticity();
                };
            };
            
            //Printing results
            printResults(iTimeStep);
        };
    };
        
    return 0;

};



//------------------------------------------------------------------------------
//----------------COMPUTE ARLEQUIN COUPLED NAVIER-STOKES PROBLEM----------------
//------------------------------------------------------------------------------
template<>
int Arlequin<2>::solveFSIArlequin(int iterNumber, double tolerance,
                                  int problem_type){

    Mat               A,F;
    Vec               b, u, All;
    PetscErrorCode    ierr;
    PetscInt          Istart, Iend, Ii, Ione, iterations;
    KSP               ksp; 
    PC                pc;
    VecScatter        ctx;
    PetscScalar       val;
    //MatNullSpace      nullsp;

    std::ofstream dragLift;
    dragLift.open("dragLift.dat", std::ofstream::out | std::ofstream::app);

    if ((problem_type < 1) || (problem_type > 2)){
        std::cout << "WRONG PROBLEM TYPE." << std::endl;
        return 0;
    };

    // Computes the system size
    int sysSize = 3 * numNodesCoarse +  
        + 3 * numNodesFine + 2 * numNodesGlueZoneFine;
    
    
    setSignaledDistance();
    setWeightFunction(1.);
    setNodalCorrespondenceFine();
    
    //STARTS NEWTON-RAPHSON
    for (int inewton = 0; inewton < iterNumber; inewton++){
        
        boost::posix_time::ptime t1 =
            boost::posix_time::microsec_clock::local_time();
        
        // Preallocates the matrix
        ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                            sysSize, sysSize, 120, NULL, 500, NULL, &A); 
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
                std::pair<Elements::LocalMatrix, Elements::LocalVector> localMV;
                localMV = elementsCoarse_[jel] -> getTransientNavierStokes();
    
                typename Elements::Connectivity connec;
                connec = elementsCoarse_[jel] -> getConnectivity();
                
                //Disperse local contributions into the global matrix
                for (int i=0; i<6; i++){
                    for (int j=0; j<6; j++){
                        
                        if (fabs(localMV.first(2*i  ,2*j  )) >= 1.e-15){
                            int dof_i = 2 * connec(i);
                            int dof_j = 2 * connec(j);
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&localMV.first(2*i  ,2*j  ),ADD_VALUES);
                        };
                        if (fabs(localMV.first(2*i+1,2*j  )) >= 1.e-15){
                            int dof_i = 2 * connec(i) + 1;
                            int dof_j = 2 * connec(j);
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&localMV.first(2*i+1,2*j  ),ADD_VALUES);
                        };

                        if (fabs(localMV.first(2*i  ,2*j+1)) >= 1.e-15){
                            int dof_i = 2 * connec(i);
                            int dof_j = 2 * connec(j) + 1;
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&localMV.first(2*i  ,2*j+1),ADD_VALUES);
                        };

                        if (fabs(localMV.first(2*i+1,2*j+1)) >= 1.e-15){
                            int dof_i = 2 * connec(i) + 1;
                            int dof_j = 2 * connec(j) + 1;
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&localMV.first(2*i+1,2*j+1),ADD_VALUES);
                        };

                        //Matrix Q and Qt
                        if (fabs(localMV.first(2*i  ,12+j)) >= 1.e-15){
                            int dof_i = 2 * connec(i);
                            int dof_j = 2 * numNodesCoarse + connec(j);
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&localMV.first(2*i  ,12+j),ADD_VALUES);
                        };
                        
                        if (fabs(localMV.first(12+j,2*i  )) >= 1.e-15){
                            int dof_i = 2 * connec(i);
                            int dof_j = 2 * numNodesCoarse + connec(j);
                            ierr = MatSetValues(A,1,&dof_j,1,&dof_i,&localMV.first(12+j,2*i  ),ADD_VALUES);
                        };

                        if (fabs(localMV.first(2*i+1,12+j)) >= 1.e-15){
                            int dof_i = 2 * connec(i) + 1;
                            int dof_j = 2 * numNodesCoarse + connec(j);
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&localMV.first(2*i+1,12+j),ADD_VALUES);
                        };

                        if (fabs(localMV.first(12+j,2*i+1)) >= 1.e-15){
                            int dof_i = 2 * connec(i) + 1;
                            int dof_j = 2 * numNodesCoarse + connec(j);
                            ierr = MatSetValues(A,1,&dof_j,1,&dof_i,&localMV.first(12+j,2*i+1),ADD_VALUES);
                        };

                        if (fabs(localMV.first(12+i,12+j)) >= 1.e-15){
                            int dof_i = 2 * numNodesCoarse + connec(i);
                            int dof_j = 2 * numNodesCoarse + connec(j);
                            ierr = MatSetValues(A,1,&dof_j,1,&dof_i,&localMV.first(12+i,12+j),ADD_VALUES);
                        };
                    };                  
                    //Rhs vector
                    if (fabs(localMV.second(2*i  )) >= 1.e-15){
                        int dof_i = 2 * connec(i);
                        ierr = VecSetValues(b,1,&dof_i,&localMV.second(2*i  ),ADD_VALUES);
                    };

                    if (fabs(localMV.second(2*i+1)) >= 1.e-15){
                        int dof_i = 2 * connec(i) + 1;
                        ierr = VecSetValues(b,1,&dof_i,&localMV.second(2*i+1),ADD_VALUES);
                    };

                    if (fabs(localMV.second(12+i)) >= 1.e-15){
                        int dof_i = 2 * numNodesCoarse + connec(i);
                        ierr = VecSetValues(b,1,&dof_i,&localMV.second(12+i),ADD_VALUES);
                    };
                };
            };
        };

        //Fine mesh
        for (int jel = 0; jel < numElemFine; jel++){   
            
            if (domDecompFine.first[jel] == rank) {
                
                //Compute Element matrix
                std::pair<Elements::LocalMatrix, Elements::LocalVector> localMV;
                localMV = elementsFine_[jel] -> getTransientNavierStokes();
    
                typename Elements::Connectivity connec;
                connec = elementsFine_[jel] -> getConnectivity();          

                //Disperse local contributions into the global matrix
                for (int i=0; i<6; i++){
                    for (int j=0; j<6; j++){
                        if (fabs(localMV.first(2*i  ,2*j  )) >= 1.e-15){
                            int dof_i = 3*numNodesCoarse + 2 * connec(i);
                            int dof_j = 3*numNodesCoarse + 2 * connec(j);
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&localMV.first(2*i  ,2*j  ),ADD_VALUES);
                        };
                        
                        if (fabs(localMV.first(2*i+1,2*j  )) >= 1.e-15){
                            int dof_i = 3*numNodesCoarse + 2 * connec(i) +1;
                            int dof_j = 3*numNodesCoarse + 2 * connec(j);
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&localMV.first(2*i+1,2*j  ),ADD_VALUES);
                        };

                        if (fabs(localMV.first(2*i  ,2*j+1)) >= 1.e-15){
                            int dof_i = 3*numNodesCoarse + 2 * connec(i);
                            int dof_j = 3*numNodesCoarse + 2 * connec(j) +1;
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&localMV.first(2*i  ,2*j+1),ADD_VALUES);
                        };

                        if (fabs(localMV.first(2*i+1,2*j+1)) >= 1.e-15){
                            int dof_i = 3*numNodesCoarse + 2 * connec(i) +1;
                            int dof_j = 3*numNodesCoarse + 2 * connec(j) +1;
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&localMV.first(2*i+1,2*j+1),ADD_VALUES);
                        };

                        //Matrix Q and Qt
                        if (fabs(localMV.first(2*i  ,12+j)) >= 1.e-15){
                            int dof_i = 3 * numNodesCoarse + 2 * connec(i);
                            int dof_j = 3 * numNodesCoarse + 2 * numNodesFine + connec(j);
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&localMV.first(2*i  ,12+j),ADD_VALUES);
                        };

                        if (fabs(localMV.first(12+j,2*i  )) >= 1.e-15){
                            int dof_i = 3 * numNodesCoarse + 2 * connec(i);
                            int dof_j = 3 * numNodesCoarse + 2 * numNodesFine + connec(j);
                            ierr = MatSetValues(A,1,&dof_j,1,&dof_i,&localMV.first(12+j,2*i  ),ADD_VALUES);
                        };

                        if (fabs(localMV.first(2*i+1,12+j)) >= 1.e-15){
                            int dof_i = 3 * numNodesCoarse + 2 * connec(i) + 1;
                            int dof_j = 3 * numNodesCoarse + 2 * numNodesFine + connec(j);
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&localMV.first(2*i+1,12+j),ADD_VALUES);
                        };
                        
                        if (fabs(localMV.first(12+j,2*i+1)) >= 1.e-15){
                            int dof_i = 3 * numNodesCoarse + 2 * connec(i) + 1;
                            int dof_j = 3 * numNodesCoarse + 2 * numNodesFine + connec(j);
                            ierr = MatSetValues(A,1,&dof_j,1,&dof_i,&localMV.first(12+j,2*i+1),ADD_VALUES);
                        };
                        
                        if (fabs(localMV.first(12+i,12+j)) >= 1.e-15){
                            int dof_i = 3 * numNodesCoarse + 2 * numNodesFine + connec(i);
                            int dof_j = 3 * numNodesCoarse + 2 * numNodesFine + connec(j);
                            ierr = MatSetValues(A,1,&dof_j,1,&dof_i,&localMV.first(12+i,12+j),ADD_VALUES);
                        };
                    };     
                    ///Rhs vector
                    if (fabs(localMV.second(2*i  )) >= 1.e-15){
                        int dof_i = 3 * numNodesCoarse + 2 * connec(i);
                        ierr = VecSetValues(b,1,&dof_i,&localMV.second(2*i  ),ADD_VALUES);
                    };

                    if (fabs(localMV.second(2*i+1)) >= 1.e-15){
                        int dof_i = 3 * numNodesCoarse + 2 * connec(i) + 1;
                        ierr = VecSetValues(b,1,&dof_i,&localMV.second(2*i+1),ADD_VALUES);
                    };

                    if (fabs(localMV.second(12+i)) >= 1.e-15){
                        int dof_i = 3 * numNodesCoarse + 2 * numNodesFine + connec(i);
                        ierr = VecSetValues(b,1,&dof_i,&localMV.second(12+i),ADD_VALUES);
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
                Rhs = - elementsFine_[jel] -> getRhsVectorLagMult();
                AStab = elementsFine_[jel] -> getJacNRMatrix();
                RhsStab = elementsFine_[jel] -> getRhsVector();
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
                            ierr = MatSetValues(A,1,&d_j,1,&d_i,
                                                &AStab(2*i  ,2*j  ),
                                                ADD_VALUES);
                        };
                        
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
                            ierr = MatSetValues(A,1,&d_j,1,&d_i,
                                                &AStab(2*i+1,2*j+1),
                                                ADD_VALUES);
                        };
                        
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
                        
                        if (fabs(AStab(2*i  ,12+j)) >= 1.e-15){
                            int dof_i = 3 * numNodesCoarse + 2 * connec(i);
                            int dof_j = 3 * numNodesCoarse + 
                                2 * numNodesFine + connec(j);
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,
                                                &AStab(2*i  ,12+j),
                                                ADD_VALUES);
                        };
                        if (fabs(AStab(2*i+1,12+j)) >= 1.e-15){
                            int dof_i = 3 * numNodesCoarse + 
                                2 * connec(i) + 1;
                            int dof_j = 3 * numNodesCoarse + 
                                2 * numNodesFine + connec(j);
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,
                                                &AStab(2*i+1,12+j),
                                                ADD_VALUES);
                        };
                    };
                    
                    //Rhs vector

                    int dof_i = 3 * numNodesCoarse + 3 * numNodesFine +
                        2 * connecL(i);
                    ierr = VecSetValues(b,1,&dof_i,&Rhs(2*i  ),
                                        ADD_VALUES);
                    
                    dof_i = 3 * numNodesCoarse + 3 * numNodesFine + 
                        2 * connecL(i) + 1;
                    ierr = VecSetValues(b,1,&dof_i,&Rhs(2*i+1),
                                        ADD_VALUES);
                    
                    dof_i = 3 * numNodesCoarse + 2 * connec(i);
                    ierr = VecSetValues(b,1,&dof_i,&rhsLagMult(2*i  )
                                        ,ADD_VALUES);
                    ierr = VecSetValues(b,1,&dof_i,&RhsStab(2*i  )
                                        ,ADD_VALUES);

                    dof_i = 3 * numNodesCoarse + 2 * connec(i) + 1;
                    ierr = VecSetValues(b,1,&dof_i,&rhsLagMult(2*i+1)
                                        ,ADD_VALUES);
                    ierr = VecSetValues(b,1,&dof_i,&RhsStab(2*i+1)
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
                    // double pspg = elementsCoarse_[iElemCoarse] -> getPSPG();
                    
                    // elementsFine_[jel] -> 
                    //     getLagrangeMultipliersDifferentMesh(iElemCoarse,pspg);
                    
                    //Computes element matrix
                    Ajac = elementsFine_[jel] -> getLagrMultMatrix();
                    Ajac = trans(Ajac);
                    AStab = elementsFine_[jel] -> getJacNRMatrix();
                    AStab = trans(AStab);
                    RhsStab = elementsFine_[jel] -> getRhsVector();
                    
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
                                ierr = MatSetValues(A,1,&dof_j,1,&dof_i,
                                                    &AStab(2*i  ,2*j  )
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
                                ierr = MatSetValues(A,1,&dof_j,1,&dof_i,
                                                    &AStab(2*i+1,2*j+1)
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
                        
                        d_i = 3 * numNodesCoarse 
                            + 3 * numNodesFine + 2 * connecL(i) + 1;
                        ierr = VecSetValues(b,1,&d_i,
                                            &lagMult.first(2*i+1),
                                            ADD_VALUES);
                        
                        int dof_i = 2 * connecC(i);
                        ierr = VecSetValues(b,1,&dof_i,
                                            &rhsLagMult(2*i  ),
                                            ADD_VALUES);
                        ierr = VecSetValues(b,1,&dof_i,
                                            &RhsStab(2*i  ),
                                            ADD_VALUES);
                        
                        dof_i = 2 * connecC(i) + 1;
                        ierr = VecSetValues(b,1,&dof_i,
                                            &rhsLagMult(2*i+1),
                                            ADD_VALUES);
                        ierr = VecSetValues(b,1,&dof_i,
                                            &RhsStab(2*i+1),
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

       
        ierr = PCFactorSetMatSolverType(pc,MATSOLVERMUMPS);
        PCFactorSetUpMatSolverType(pc);
        PCFactorGetMatrix(pc,&F);

        PetscInt ival,icntl;
        icntl = 14; ival = 60;
        MatMumpsSetIcntl(F,icntl,ival);
        
#endif
        
        
        ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
        ierr = KSPSetUp(ksp);
        
        
#if defined(PETSC_HAVE_MUMPS)
        PetscInt  info1,info2,icntl14;

        MatMumpsGetInfo(F,1,&info1);
        MatMumpsGetInfo(F,2,&info2);
        MatMumpsGetIcntl(F,14,&icntl14);
        if(rank == 0) std::cout << " INFO(1) = " << info1
                                << " " << info2 << " " << icntl14 << std::endl;
#endif


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
        nodesFine_[i] -> setVelocityArlequin(0,nodesFine_[i] -> getVelocity(0));
        nodesFine_[i] -> setVelocityArlequin(1,nodesFine_[i] -> getVelocity(1));
        nodesFine_[i] -> setPressureArlequin(nodesFine_[i] -> getPressure());
    };
    
    for (int i = 0; i<numNodesGlueZoneFine; i++){
        typename Nodes::VecLocD xsi;
        typename Quadrature::NodalValuesQuad u_coarse, v_coarse, p_coarse;
        typename Elements::Connectivity connecCoarse;
        
        double u = 0.;
        double v = 0.;
        double p = 0.;
        
        int elCoarse = nodesFine_[nodesGlueZoneFine_[i]] -> 
            getNodalElemCorrespondence();
        xsi = nodesFine_[nodesGlueZoneFine_[i]] -> getNodalXsiCorrespondence();
        
        connecCoarse = elementsCoarse_[elCoarse] -> getConnectivity();
        
        for (int j=0; j<6; j++){
            u_coarse(j) = nodesCoarse_[connecCoarse(j)] -> getVelocity(0);
            v_coarse(j) = nodesCoarse_[connecCoarse(j)] -> getVelocity(1);
            p_coarse(j) = nodesCoarse_[connecCoarse(j)] -> getPressure();
        };
        
        shapeQuad.evaluate(xsi,phi_);
        
        for (int j=0; j<6; j++){
            u += u_coarse(j) * phi_(j);
            v += v_coarse(j) * phi_(j);
            p += p_coarse(j) * phi_(j);
        };

        
        double wFunc = nodesFine_[nodesGlueZoneFine_[i]] -> getWeightFunction();
        
        double u_int = nodesFine_[nodesGlueZoneFine_[i]] -> 
            getVelocity(0) * wFunc + u * (1. - wFunc);
        double v_int = nodesFine_[nodesGlueZoneFine_[i]] -> 
            getVelocity(1) * wFunc + v * (1. - wFunc);
        double p_int = nodesFine_[nodesGlueZoneFine_[i]] -> 
            getPressure() * wFunc + p * (1. - wFunc);
        
        //std::cout << "asdasd " << elCoarse << " " << p_int << " " << nodesFine_[nodesGlueZoneFine_[i]] -> getPressure() << " " << wFunc <<  std::endl;


        nodesFine_[nodesGlueZoneFine_[i]] -> setVelocityArlequin(0,u_int);
        nodesFine_[nodesGlueZoneFine_[i]] -> setVelocityArlequin(1,v_int);
        nodesFine_[nodesGlueZoneFine_[i]] -> setPressureArlequin(p_int);
        
    };
    
    for (int i=0; i<numNodesCoarse; i++){
        nodesCoarse_[i] -> 
            setVelocityArlequin(0,nodesCoarse_[i] -> getVelocity(0));
        nodesCoarse_[i] -> 
            setVelocityArlequin(1,nodesCoarse_[i] -> getVelocity(1));
        nodesCoarse_[i] -> 
            setPressureArlequin(nodesCoarse_[i] -> getPressure());
    };
      
    
    return 0;

};


#endif
