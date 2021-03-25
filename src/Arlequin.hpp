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

template<int DIM, int DEG>
class Arlequin{
public:
    /// Defines the class Fluid locally
    typedef Fluid<DIM,DEG>                 FluidMesh;

    /// Defines the class Element locally
    typedef typename FluidMesh::Elements   Elements;

    /// Defines the class Node locally
    typedef typename FluidMesh::Node       Nodes;

    /// Defines the class Boundary locally
    typedef typename FluidMesh::Boundaries Boundary;

    /// Defines the class SpecialQuad locally
    typedef typename Elements::SpecialQuad Quadrature;

    /// Defines the class Glue locally
    typedef Glue<DIM,DEG>                  GlueZone;

    typedef FluidParameters<DIM,DEG>       Parameters;

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

    Mat               A, F;
    Vec               b, u, All;
    PetscErrorCode    ierr;
    PetscInt          Istart, Iend, Ii, Ione, iterations;
    KSP               ksp; 
    PC                pc;
    VecScatter        ctx;
    PetscScalar       val;
    PetscLogDouble bytes = 0;

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
                         int problem_type, int iTimeStep);

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
    void searchNodeCorrespondence(double* x, std::vector<Nodes *> nodes,
                                  std::vector<Elements *> elements,
                                  int numElem, int &elCorr, double* xsiCorr, int elSearch);

    void setMatVecValuesFineModel(double **matrix, double *rhs, int* connec);
    void setMatVecValuesCoarseModel(double **matrix, double *rhs, int* connec);
    void setMatVecValuesLagMultFineFine(double **Ajac2, double **localMV_mat, 
                                        double **ArlequinA1, double **ArlequinA2, 
                                        double *Rhs2, double *rhsLagMult2,
                                        double *localMV_vec, double *RhsArlequin2,
                                        int* connec, int* connecL);

    void setMatVecValuesLagMultFineCoarse(double **Ajac2, double **localMV_mat, 
                                          double **ArlequinA1, double **ArlequinA2, 
                                          double *Rhs2, double *rhsLagMult2,
                                          double *localMV_vec, double *RhsArlequin2,
                                          int* connecC, int* connecL);

    void assembleArlequinSystem();

    /// Print the results for Paraview post-processing
    /// @param int time step
    void printResultsCoarse(int step);
    void printResultsFine(int step);

    void initialAcceleration();

};

//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//-------------------------COMPUTE ELEMENT REGIONS/BOXES------------------------
//------------------------------------------------------------------------------
template<>
void Arlequin<2,2>::setElementBoxes() {
    
    int *connec;
    double xk[2], Xk[2];

    //Compute element boxes for coarse model
    //Only function for straight elements
    for (int jel = 0; jel < numElemCoarse; jel++){
        connec = elementsCoarse_[jel] -> getConnectivity();
        double* x1 = nodesCoarse_[connec[0]] -> getCoordinates();
        double* x2 = nodesCoarse_[connec[1]] -> getCoordinates();
        double* x3 = nodesCoarse_[connec[2]] -> getCoordinates();      

        xk[0] = std::min(x1[0],std::min(x2[0], x3[0]));
        xk[1] = std::min(x1[1],std::min(x2[1], x3[1]));

        Xk[0] = std::max(x1[0],std::max(x2[0], x3[0]));
        Xk[1] = std::max(x1[1],std::max(x2[1], x3[1]));        
        
        elementsCoarse_[jel] -> setIntersectionParameters(xk, Xk);
    };

    //Compute element boxes for fine model
    //Only function for straight elements
    for (int jel = 0; jel < numElemFine; jel++){
        connec = elementsFine_[jel] -> getConnectivity();
        double* x1 = nodesFine_[connec[0]] -> getCoordinates();
        double* x2 = nodesFine_[connec[1]] -> getCoordinates();
        double* x3 = nodesFine_[connec[2]] -> getCoordinates();      

        xk[0] = std::min(x1[0],std::min(x2[0], x3[0]));
        xk[1] = std::min(x1[1],std::min(x2[1], x3[1]));

        Xk[0] = std::max(x1[0],std::max(x2[0], x3[0]));
        Xk[1] = std::max(x1[1],std::max(x2[1], x3[1]));        
        
        elementsFine_[jel] -> setIntersectionParameters(xk, Xk);
    };

    return;
};


//------------------------------------------------------------------------------
//-------------------COMPUTE NODAL CORRESPONDECE WITH ELEMENTS------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Arlequin<DIM,DEG>::searchNodeCorrespondence(double* x,std::vector<Nodes *> nodes, 
                                           std::vector<Elements *> elements, 
                                           int numElem, int &elCorr, double* xsiCorr, int elSearch){
    
    
    int *connec;
    QuadShapeFunction<2,2>                       shapeQuad;
    double phi_[6] = {};

    double **ainv;
    ainv = new double*[DIM];
    for (int i = DIM; i--; ) ainv[i] = new double[DIM];
    
    double xsiCC[3];
    std::pair<double*,double*> XK;
    std::pair<double*,double*> dCk;

    elCorr = 150000;
    xsiCorr[0] = 1.e50;
    xsiCorr[1] = 1.e50;
    double xsi[2];
    double x_[2];
    double deltaX[2] = {};
    double deltaXsi[2] = {};
    bool flg = true;

    // std::cout << "elSearch " << elSearch << std::endl;


    connec = elements[elSearch] -> getConnectivity();
    xsiCC[0] = 1.e10;
    xsiCC[1] = 1.e10;
    xsiCC[2] = 1.e10;
    
    xsi[0] = 1. / 3.;
    xsi[1] = 1. / 3.;

    shapeQuad.evaluate(xsi,phi_);
    
    x_[0] = 0.;
    x_[1] = 0.;

    for (int i = 0; i < 6; i++){
        double* xint = nodes[connec[i]] -> getCoordinates();
        x_[0] += xint[0] * phi_[i];
        x_[1] += xint[1] * phi_[i];                    
    };

    double error = 1.e6;
    
    int iterations = 0;

    while ((error > 1.e-8) && (iterations < 4)) {
        
        iterations++;
        
        deltaX[0] = x[0] - x_[0];
        deltaX[1] = x[1] - x_[1];
        
        deltaXsi[0] = 0.;
        deltaXsi[1] = 0.;
        
        elements[elSearch] -> getJacobianMatrix(xsi,ainv);
    
        // noalias(deltaXsi) = prod(trans(ainv),deltaX);

        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
                deltaXsi[i] += ainv[j][i] * deltaX[j];
        
        xsi[0] += deltaXsi[0];
        xsi[1] += deltaXsi[1];
        
        x_[0] = 0.;
        x_[1] = 0.;
    
        shapeQuad.evaluate(xsi,phi_);
        
        for (int i=0; i<6; i++){
            double* xint = nodes[connec[i]] -> getCoordinates();
            x_[0] += xint[0] * phi_[i];
            x_[1] += xint[1] * phi_[i];                    
        };                   
        error = sqrt(deltaXsi[0]*deltaXsi[0] + deltaXsi[1]*deltaXsi[1]);
    };
    
    double t1 = -1.e-2;
    double t2 =  1. - t1;
    
    xsiCC[0] = xsi[0];
    xsiCC[1] = xsi[1];       
    xsiCC[2] = 1. - xsiCC[0] - xsiCC[1];

    if ((xsiCC[0] >= t1) && (xsiCC[1] >= t1) && (xsiCC[2] >= t1) &&
        (xsiCC[0] <= t2) && (xsiCC[1] <= t2) && (xsiCC[2] <= t2)){

        xsiCorr[0] = xsi[0]; xsiCorr[1] = xsi[1];
        elCorr = elSearch;
        // return;
    } else {
        for (int jel = 0; jel < numElem; jel++){
            connec = elements[jel] -> getConnectivity();

            //get boxes information        
            XK = elements[jel] -> getXIntersectionParameter();

            //Chech if the node is inside the element box
            if ((x[0] < XK.first[0]) || (x[0] > XK.second[0]) ||
                (x[1] < XK.first[1]) || (x[1] > XK.second[1])) continue;
            
            //Compute nodal correspondence
            xsiCC[0] = 1.e10;
            xsiCC[1] = 1.e10;
            xsiCC[2] = 1.e10;
            
            xsi[0] = 1. / 3.;
            xsi[1] = 1. / 3.;

            shapeQuad.evaluate(xsi,phi_);
            
            x_[0] = 0.;
            x_[1] = 0.;

            for (int i = 0; i < 6; i++){
                double* xint = nodes[connec[i]] -> getCoordinates();
                x_[0] += xint[0] * phi_[i];
                x_[1] += xint[1] * phi_[i];                    
            };

            double error = 1.e6;
            
            int iterations = 0;

            while ((error > 1.e-8) && (iterations < 4)) {
                
                iterations++;
                
                deltaX[0] = x[0] - x_[0];
                deltaX[1] = x[1] - x_[1];
                
                deltaXsi[0] = 0.;
                deltaXsi[1] = 0.;
                
                elements[jel] -> getJacobianMatrix(xsi,ainv);
            
                // noalias(deltaXsi) = prod(trans(ainv),deltaX);

                for (int i = 0; i < 2; i++)
                    for (int j = 0; j < 2; j++)
                        deltaXsi[i] += ainv[j][i] * deltaX[j];
                
                xsi[0] += deltaXsi[0];
                xsi[1] += deltaXsi[1];
                
                x_[0] = 0.;
                x_[1] = 0.;
            
                shapeQuad.evaluate(xsi,phi_);
                
                for (int i=0; i<6; i++){
                    double* xint = nodes[connec[i]] -> getCoordinates();
                    x_[0] += xint[0] * phi_[i];
                    x_[1] += xint[1] * phi_[i];                    
                };                   
                error = sqrt(deltaXsi[0]*deltaXsi[0] + deltaXsi[1]*deltaXsi[1]);
            };
            
            double t1 = -1.e-2;
            double t2 =  1. - t1;
            
            xsiCC[0] = xsi[0];
            xsiCC[1] = xsi[1];       
            xsiCC[2] = 1. - xsiCC[0] - xsiCC[1];

            if ((xsiCC[0] >= t1) && (xsiCC[1] >= t1) && (xsiCC[2] >= t1) &&
                (xsiCC[0] <= t2) && (xsiCC[1] <= t2) && (xsiCC[2] <= t2)){

                xsiCorr[0] = xsi[0]; xsiCorr[1] = xsi[1];
                elCorr = jel;
                break;
            }
        }
    };

    for (int i = DIM; i--; ) delete [] ainv[i];
    delete [] ainv;

    if (fabs(xsi[0]) > 2.) std::cout << "PROBEM SEARCHING NODE CORRESPONDENCE " 
                                     << std::endl;

    return;
};

//------------------------------------------------------------------------------
//--------COMPUTE NODAL CORRESPONDECE OF FINE NODES WITH COARSE ELEMENTS--------
//------------------------------------------------------------------------------
template<>
void Arlequin<2,2>::setNodalCorrespondenceFine() {

    //FINE MESH
    // for (int inode=0; inode < numNodesFine; inode ++){
    //     nodesFine_[inode] -> setNodalCorrespondence(1.e10,xsi);
    // };

    double &alpha_f = parametersFine -> getAlphaF();

    for (int inode = 0; inode < numNodesGlueZoneFine; inode++) {
        
        double* x = nodesFine_[nodesGlueZoneFine_[inode]] -> getCoordinates();

        int elCorr = 0;
        double xsiCorr[2] = {};


        // searchNodeCorrespondence(x, nodesCoarse_, elementsCoarse_[nodesFine_[nodesGlueZoneFine_[inode]]->getNodalElemCorrespondence()], 
        //                          1,elCorr,xsiCorr, );

        // double t1 = -1.e-2;
        // double t2 =  1. - t1;

        // if ((xsiCorr[0] >= t1) && (xsiCorr[1] >= t1) && (1. - xsiCorr[0] - xsiCorr[1] >= t1) &&
        //     (xsiCorr[0] <= t2) && (xsiCorr[1] <= t2) && (1. - xsiCorr[0] - xsiCorr[1] <= t2)){
            
        //     nodesFine_[nodesGlueZoneFine_[inode]] -> setNodalCorrespondence(elCorr,xsiCorr);
        // } else {
            searchNodeCorrespondence(x, nodesCoarse_, elementsCoarse_, 
                                 elementsCoarse_.size(),elCorr,xsiCorr,nodesFine_[nodesGlueZoneFine_[inode]] -> getNodalElemCorrespondence());
            nodesFine_[nodesGlueZoneFine_[inode]] -> setNodalCorrespondence(elCorr,xsiCorr);             
        // }

        
   
        // std::cout << "CORRESP " << elCorr << " " << corresp.first << std::endl 
                  // << xsiCorr[0] << " " << xsiCorr[1] << " " << corresp.second[0] << " " << corresp.second[1] << std::endl;

        
            
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
        
        double x1[6], x2[6];
        int *connec;

        connec = elementsFine_[elementsGlueZoneFine_[ielem]] -> getConnectivity();
        
        for (int i = 0; i < 6; i++){
            double *x = nodesFine_[connec[i]] -> getCoordinates();
            double *xp = nodesFine_[connec[i]] -> getPreviousCoordinates();
            
            x1[i] = alpha_f * x[0] + (1. - alpha_f) * xp[0];
            x2[i] = alpha_f * x[1] + (1. - alpha_f) * xp[1];
        };

        // std::cout << "XX1 " << x2 << " " << x22 << " " << x222 << std::endl;

        int numberIntPoints = elementsFine_[elementsGlueZoneFine_[ielem]] -> 
            getNumberOfIntegrationPoints();

        for (int i = 0; i < numberIntPoints; i++){

            // if (elementsFine_[ielem] -> getIntegPointInGlueZone(i) == true){
                double x[2];

                x[0] = quad.interpolateQuadraticVariable(x1,i);
                x[1] = quad.interpolateQuadraticVariable(x2,i);

                int elCorr = 0;
                double xsiCorr[2] = {};
            
                searchNodeCorrespondence(x,nodesCoarse_,elementsCoarse_,
                                         elementsCoarse_.size(),elCorr,xsiCorr,elementsFine_[elementsGlueZoneFine_[ielem]] -> getIntegPointCorrespondenceElement(i));
                // if(ielem == 0){
                //     if(i==0)std::cout <<" COOr " << corresp.first << " " << corresp.second(0) << " " << corresp.second(1) << std::endl;
                // };

                 // std::cout << "CORRESP " << elCorr << " " << corresp.first << std::endl 
                 //  << xsiCorr[0] << " " << xsiCorr[1] << " " << corresp.second[0] << " " << corresp.second[1] << std::endl;

                elementsFine_[elementsGlueZoneFine_[ielem]] -> 
                    setIntegrationPointCorrespondence(i,elCorr,xsiCorr);
                // };
        };

    };
};

//------------------------------------------------------------------------------
//--------COMPUTE NODAL CORRESPONDECE OF COARSE NODES WITH FINE ELEMENTS--------
//------------------------------------------------------------------------------
template<>
void Arlequin<2,2>::setNodalCorrespondenceCoarse() {

    //COARSE MESH
    for (int inode=0; inode < numNodesCoarse; inode ++){
        double xsi[2] = {};
        nodesCoarse_[inode] -> setNodalCorrespondence(1.e10,xsi);
    };

    for (int inode = 0; inode < numNodesGlueZoneFine; inode++) {
        
        double *x = nodesCoarse_[nodesGlueZoneFine_[inode]] -> getCoordinates();
        
        int elCorr = 0;
        double xsiCorr[2] = {};

        searchNodeCorrespondence(x,nodesFine_,elementsFine_, 
                                 elementsFine_.size(),elCorr,xsiCorr,0);
   
        nodesCoarse_[nodesGlueZoneFine_[inode]] -> setNodalCorrespondence(elCorr,xsiCorr);
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
        
        double x1[6], x2[6];
        int *connec;

        connec = elementsCoarse_[elementsGlueZoneFine_[ielem]] -> getConnectivity();
        
        for (int i = 0; i < 6; i++){
            double *x = nodesCoarse_[connec[i]] -> getCoordinates();
            
            x1[i] = x[0];
            x2[i] = x[1];
        };

        int numberIntPoints = elementsCoarse_[elementsGlueZoneFine_[ielem]] -> 
            getNumberOfIntegrationPoints();

        for (int i = 0; i < numberIntPoints; i++){

            if (elementsCoarse_[ielem] -> getIntegPointInGlueZone(i) == true){
                
                double x[2];

                x[0] = quad.interpolateQuadraticVariable(x1,i);
                x[1] = quad.interpolateQuadraticVariable(x2,i);
                
                int elCorr = 0;
                double xsiCorr[2] = {};

                searchNodeCorrespondence(x,nodesFine_,elementsFine_,
                                         elementsFine_.size(),elCorr,xsiCorr,0);

                elementsCoarse_[elementsGlueZoneFine_[ielem]] -> 
                    setIntegrationPointCorrespondence(i,elCorr,xsiCorr);
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
void Arlequin<2,2>::setSignaledDistance(){


    int *connec;
    int bconnec[3];
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
        if (boundaryFine_[i]->getConstrain(0) >= 2){
            // std::cout <<"AQUIISDUAS  "  << rank << std::endl;
            connec = elementsFine_[boundaryFine_[i] -> getElement()] -> getConnectivity();
                            
            if (elementsFine_[boundaryFine_[i] -> getElement()] -> getElemSideInBoundary() == 0){
                bconnec[0] = connec[1];
                bconnec[1] = connec[2];
                bconnec[2] = connec[4];
            };
            if (elementsFine_[boundaryFine_[i] -> getElement()] -> getElemSideInBoundary() == 1){
                bconnec[0] = connec[2];
                bconnec[1] = connec[0];
                bconnec[2] = connec[5];
            };
            if (elementsFine_[boundaryFine_[i] -> getElement()] -> getElemSideInBoundary() == 2){
                bconnec[0] = connec[0];
                bconnec[1] = connec[1];
                bconnec[2] = connec[3];
            };
            
            //bconnec = boundaryFine_[i]->getBoundaryConnectivity();
            
            //first segment
            int no1 = bconnec[0];
            int no2 = bconnec[2];
            double *x1 = nodesFine_[no1] -> getCoordinates();
            double *x2 = nodesFine_[no2] -> getCoordinates();

            double sLength = sqrt((x2[1] - x1[1]) * (x2[1] - x1[1]) +
                                  (x1[0] - x2[0]) * (x1[0] - x2[0]));

            double n[2];
            n[0] = (x2[1] - x1[1]) / sLength;
            n[1] = (x1[0] - x2[0]) / sLength;

            nodesFine_[no1] -> setInnerNormal(n);
            nodesFine_[no2] -> setInnerNormal(n);

            //second segment
            no1 = bconnec[2];
            no2 = bconnec[1];
            x1 = nodesFine_[no1] -> getCoordinates();
            x2 = nodesFine_[no2] -> getCoordinates();

            sLength = sqrt((x2[1] - x1[1]) * (x2[1] - x1[1]) +
                           (x1[0] - x2[0]) * (x1[0] - x2[0]));

            n[0] = (x2[1] - x1[1]) / sLength;
            n[1] = (x1[0] - x2[0]) / sLength;

            nodesFine_[no1] -> setInnerNormal(n);
            nodesFine_[no2] -> setInnerNormal(n);
            
        };
    };


    //Fine mesh nodes
    for (int ino = 0; ino < numNodesFine; ino++){
        double* x = nodesFine_[ino] -> getCoordinates();
        dist=10000000000000000000000000000.;
        
        for (int i = 0; i < numBoundElemFine; i++){
            if (boundaryFine_[i] -> getConstrain(0) == 2){
                connec = elementsFine_[boundaryFine_[i] -> getElement()] -> getConnectivity();
                
                if (elementsFine_[boundaryFine_[i] -> getElement()] -> getElemSideInBoundary() == 0){
                    bconnec[0] = connec[1];
                    bconnec[1] = connec[2];
                    bconnec[2] = connec[4];
                };
                if (elementsFine_[boundaryFine_[i] -> getElement()] -> getElemSideInBoundary() == 1){
                    bconnec[0] = connec[2];
                    bconnec[1] = connec[0];
                    bconnec[2] = connec[5];
                };
                if (elementsFine_[boundaryFine_[i] -> getElement()] -> getElemSideInBoundary() == 2){
                    bconnec[0] = connec[0];
                    bconnec[1] = connec[1];
                    bconnec[2] = connec[3];
                };
                //bconnec = boundaryFine_[i] -> getBoundaryConnectivity();

                //first segment
                int no1 = bconnec[0];
                int no2 = bconnec[2];
                // std::cout<<no1<<" nos "<<no2<<std::endl;
              
                double* x1 = nodesFine_[no1] -> getCoordinates();
                double* x2 = nodesFine_[no2] -> getCoordinates();
                
                double aux0 =  sqrt((x2[1] - x1[1]) * (x2[1] - x1[1]) +
                                    (x2[0] - x1[0]) * (x2[0] - x1[0]));
                double aux1 = ((x[0] - x1[0]) * (x2[0] - x1[0])+
                               (x[1] - x1[1]) * (x2[1] - x1[1])) / aux0;
                double dist2 =-((x2[1] - x1[1]) * x[0] - 
                                (x2[0] - x1[0]) * x[1] +
                                x2[0] * x1[1] - x2[1] * x1[0]) / aux0;
                
                if (aux1 > aux0){
                    dist2 = sqrt((x2[1] - x[1]) * (x2[1] - x[1]) +
                                 (x2[0] - x[0]) * (x2[0] - x[0]));
                    //find signal
                    //side normal vector
                    double *n = nodesFine_[no2] -> getInnerNormal();
                    double test[2];

                    test[0] = x[0] - x2[0];
                    test[1] = x[1] - x2[1];
                    double signaltest = n[0]*test[0] + n[1]*test[1];
                    double signal = -1.;
                    
                    if (signaltest <= -0.001)signal = 1.;
                    
                    dist2 *= signal;
                };

                if (aux1 < 0.){
                    dist2 = sqrt((x[1] - x1[1]) * (x[1] - x1[1]) +
                                 (x[0] - x1[0]) * (x[0] - x1[0]));
                    //find signal
                    //side normal vector
                    double *n = nodesFine_[no1] -> getInnerNormal();
                    double test[2];

                    test[0] = x[0] - x1[0];
                    test[1] = x[1] - x1[1];
                    double signaltest = n[0]*test[0] + n[1]*test[1];
                    double signal = -1.;
                    
                    if (signaltest <= -0.001) signal = 1.;
                    
                    dist2 *= signal;
                };
                
                if (fabs(dist2) < fabs(dist)) dist = dist2;
                
                //second segment
                no1 = bconnec[2];
                no2 = bconnec[1];
                x1 = nodesFine_[no1] -> getCoordinates();
                x2 = nodesFine_[no2] -> getCoordinates();
                
                aux0 = sqrt((x2[1] - x1[1]) * (x2[1] - x1[1]) +
                            (x2[0] - x1[0]) * (x2[0] - x1[0]));
                aux1 = ((x[0] - x1[0]) * (x2[0] - x1[0]) + 
                        (x[1] - x1[1]) * (x2[1] - x1[1])) / aux0;
                dist2 = -((x2[1] - x1[1]) * x[0] - (x2[0] - x1[0]) * x[1] +
                          x2[0] * x1[1] - x2[1] * x1[0]) / aux0;

                if (aux1 > aux0){
                    dist2 = sqrt((x2[1] - x[1]) * (x2[1] - x[1]) +
                                 (x2[0] - x[0]) * (x2[0] - x[0]));
                    double *n = nodesFine_[no2] -> getInnerNormal();
                    
                    double test[2];
                    test[0] = x[0] - x2[0];
                    test[1] = x[1] - x2[1];
                    double signaltest = n[0]*test[0] + n[1]*test[1];
                    double signal = -1.;
                    
                    if (signaltest <= -0.001)signal = 1.;
                    
                    dist2 *= signal;
                                       
                };

                if (aux1 < 0.){
                    dist2 = sqrt((x[1] - x1[1]) * (x[1] - x1[1]) +
                                 (x[0] - x1[0]) * (x[0] - x1[0]));
                    //find signal
                    //side normal vector
                    double *n = nodesFine_[no1] -> getInnerNormal();
                    
                    double test[2];
                    test[0] = x[0] - x1[0];
                    test[1] = x[1] - x1[1];
                    double signaltest = n[0]*test[0] + n[1]*test[1];
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
        double* x = nodesCoarse_[ino] -> getCoordinates();
        dist=10000000000000000000000000000.;
        
        for (int i = 0; i < numBoundElemFine; i++){
            if (boundaryFine_[i] -> getConstrain(0) == 2){

                connec = elementsFine_[boundaryFine_[i] -> getElement()] -> getConnectivity();
                
                if (elementsFine_[boundaryFine_[i] -> getElement()] -> getElemSideInBoundary() == 0){
                    bconnec[0] = connec[1];
                    bconnec[1] = connec[2];
                    bconnec[2] = connec[4];
                };
                if (elementsFine_[boundaryFine_[i] -> getElement()] -> getElemSideInBoundary() == 1){
                    bconnec[0] = connec[2];
                    bconnec[1] = connec[0];
                    bconnec[2] = connec[5];
                };
                if (elementsFine_[boundaryFine_[i] -> getElement()] -> getElemSideInBoundary() == 2){
                    bconnec[0] = connec[0];
                    bconnec[1] = connec[1];
                    bconnec[2] = connec[3];
                };
                //bconnec = boundaryFine_[i] -> getBoundaryConnectivity();

                //first segment
                int no1 = bconnec[0];
                int no2 = bconnec[2];
                // std::cout<<no1<<" nos "<<no2<<std::endl;
              
                double* x1 = nodesFine_[no1] -> getCoordinates();
                double* x2 = nodesFine_[no2] -> getCoordinates();
                
                double aux0 =  sqrt((x2[1] - x1[1]) * (x2[1] - x1[1]) +
                                    (x2[0] - x1[0]) * (x2[0] - x1[0]));
                double aux1 = ((x[0] - x1[0]) * (x2[0] - x1[0])+
                               (x[1] - x1[1]) * (x2[1] - x1[1])) / aux0;
                double dist2 =-((x2[1] - x1[1]) * x[0] - 
                                (x2[0] - x1[0]) * x[1] +
                                x2[0] * x1[1] - x2[1] * x1[0]) / aux0;
                
                if (aux1 > aux0){
                    dist2 = sqrt((x2[1] - x[1]) * (x2[1] - x[1]) +
                                 (x2[0] - x[0]) * (x2[0] - x[0]));
                    //find signal
                    //side normal vector
                    double *n = nodesFine_[no2] -> getInnerNormal();
                    
                    double test[2];
                    test[0] = x[0] - x2[0];
                    test[1] = x[1] - x2[1];
                    double signaltest = n[0]*test[0] + n[1]*test[1];
                    double signal = -1.;
                    
                    if (signaltest <= -0.001)signal = 1.;
                    
                    dist2 *= signal;
                };

                if (aux1 < 0.){
                    dist2 = sqrt((x[1] - x1[1]) * (x[1] - x1[1]) +
                                 (x[0] - x1[0]) * (x[0] - x1[0]));
                    //find signal
                    //side normal vector
                    double *n = nodesFine_[no1] -> getInnerNormal();
                    
                    double test[2];
                    test[0] = x[0] - x1[0];
                    test[1] = x[1] - x1[1];
                    double signaltest = n[0]*test[0] + n[1]*test[1];
                    double signal = -1.;
                    
                    if (signaltest <= -0.001) signal = 1.;
                    
                    dist2 *= signal;
                };
                
                if (fabs(dist2) < fabs(dist)) dist = dist2;
                
                //second segment
                no1 = bconnec[2];
                no2 = bconnec[1];
                x1 = nodesFine_[no1] -> getCoordinates();
                x2 = nodesFine_[no2] -> getCoordinates();
                
                aux0 = sqrt((x2[1] - x1[1]) * (x2[1] - x1[1]) +
                            (x2[0] - x1[0]) * (x2[0] - x1[0]));
                aux1 = ((x[0] - x1[0]) * (x2[0] - x1[0]) + 
                        (x[1] - x1[1]) * (x2[1] - x1[1])) / aux0;
                dist2 = -((x2[1] - x1[1]) * x[0] - (x2[0] - x1[0]) * x[1] +
                          x2[0] * x1[1] - x2[1] * x1[0]) / aux0;

                if (aux1 > aux0){
                    dist2 = sqrt((x2[1] - x[1]) * (x2[1] - x[1]) +
                                 (x2[0] - x[0]) * (x2[0] - x[0]));
                    double *n = nodesFine_[no2] -> getInnerNormal();
                    
                    double test[2];
                    test[0] = x[0] - x2[0];
                    test[1] = x[1] - x2[1];
                    double signaltest = n[0]*test[0] + n[1]*test[1];
                    double signal = -1.;
                    
                    if (signaltest <= -0.001)signal = 1.;
                    
                    dist2 *= signal;
                                       
                };

                if (aux1 < 0.){
                    dist2 = sqrt((x[1] - x1[1]) * (x[1] - x1[1]) +  
                                 (x[0] - x1[0]) * (x[0] - x1[0]));
                    //find signal
                    //side normal vector
                    double *n = nodesFine_[no1] -> getInnerNormal();
                    
                    double test[2];
                    test[0] = x[0] - x1[0];
                    test[1] = x[1] - x1[1];
                    double signaltest = n[0]*test[0] + n[1]*test[1];
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
void Arlequin<2,2>::setCouplingZone(){

    int *connec;

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
            double* x = nodesFine_[connec[ino]] -> getCoordinates();
            double dist = nodesFine_[connec[ino]] -> getDistFunction();
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
            nodesCZ[connec[ino]] += 1;
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
        double* x = nodesFine_[nodesGlueZoneFine_[i]] -> getCoordinates();
        
        Nodes *no = new Nodes(x,i);
        nodesLagrangeFine_.push_back(no);
    };

    for (int i = 0; i < numElemGlueZoneFine; i++){
        int *connecAux;
        connecAux = new int[6];

        connec = elementsFine_[elementsGlueZoneFine_[i]] -> getConnectivity();
        
        for (int ino = 0; ino < numNodesGlueZoneFine; ino++){
            if (nodesGlueZoneFine_[ino] == connec[0]) connecAux[0] = ino;
            if (nodesGlueZoneFine_[ino] == connec[1]) connecAux[1] = ino;
            if (nodesGlueZoneFine_[ino] == connec[2]) connecAux[2] = ino;
            if (nodesGlueZoneFine_[ino] == connec[3]) connecAux[3] = ino;
            if (nodesGlueZoneFine_[ino] == connec[4]) connecAux[4] = ino;
            if (nodesGlueZoneFine_[ino] == connec[5]) connecAux[5] = ino;
        };
        
        glueZoneFine_[i] -> setConnectivity(connecAux);
        // glueZoneFine_[i] -> setNodes(nodesLagrangeFine_);

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
            double* x = nodesCoarse_[connec[ino]] -> getCoordinates();

            if ((x[0] < lim1) || (x[0] > lim2) || 
                (x[1] < lim1) || (x[1] > lim2)){

            }else{
                if ((x[0] > lim1 + tick) && (x[0] < lim2 - tick) &&
                    (x[1] > lim1 + tick) && (x[1] < lim2 - tick)){

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
            
            for (int i=0; i < elementsCoarse_[jel] -> getNumberOfIntegrationPoints(); i++){
                double* x = elementsCoarse_[jel] -> getIntegPointCoordinatesValue(i);

                if ((x[0] < lim1) || (x[0] > lim2) || 
                    (x[1] < lim1) || (x[1] > lim2)){

                }else{
                    if ((x[0] > lim1 + tick) && (x[0] < lim2 - tick) &&
                        (x[1] > lim1 + tick) && (x[1] < lim2 - tick)){

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
            nodesCZ2[connec[ino]] += 1;
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
        double* x = nodesCoarse_[nodesGlueZoneCoarse_[i]] -> getCoordinates();
        
        Nodes *no = new Nodes(x,i);
        nodesLagrangeCoarse_.push_back(no);
    };

    for (int i = 0; i < numElemGlueZoneCoarse; i++){
        int *connecAux;
        connecAux = new int[6];

        connec = elementsCoarse_[elementsGlueZoneCoarse_[i]] -> getConnectivity();
        
        for (int ino = 0; ino < numNodesGlueZoneCoarse; ino++){
            if (nodesGlueZoneCoarse_[ino] == connec[0]) connecAux[0] = ino;
            if (nodesGlueZoneCoarse_[ino] == connec[1]) connecAux[1] = ino;
            if (nodesGlueZoneCoarse_[ino] == connec[2]) connecAux[2] = ino;
            if (nodesGlueZoneCoarse_[ino] == connec[3]) connecAux[3] = ino;
            if (nodesGlueZoneCoarse_[ino] == connec[4]) connecAux[4] = ino;
            if (nodesGlueZoneCoarse_[ino] == connec[5]) connecAux[5] = ino;
        };
        
        glueZoneCoarse_[i] -> setConnectivity(connecAux);
        // glueZoneCoarse_[i] -> setNodes(nodesLagrangeCoarse_);

    };
     
};

//------------------------------------------------------------------------------
//----------------------COMPUTES THE WEIGHT FUNCTION VALUE----------------------
//------------------------------------------------------------------------------
template<>
double Arlequin<2,2>::weightFunctionFineValue(double r, double epsilon){

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
double Arlequin<2,2>::weightFunctionCoarseValue(double r, double epsilon){

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
void Arlequin<2,2>::setWeightFunction(double val){
    
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
//-----------------------------PRINT COARSE RESULTS-----------------------------
//------------------------------------------------------------------------------
template<>
void Arlequin<2,2>::printResultsCoarse(int step) {

    //PRINT COARSE MODEL RESULTS
    std::string result;
    std::ostringstream convert;

    convert << step+100000;
    result = convert.str();

    std::string s = "saidaVelCoarse"+result+".xdmf";
    std::fstream output_v(s.c_str(), std::ios_base::out);

    std::string s1 = "resultCoarse"+result+".h5";
    std::fstream filename(s1.c_str(), std::ios_base::out);

    std::string s2 = "geometryCoarse.h5";
    if (step == 0){
        std::fstream filename(s2.c_str(), std::ios_base::out);
    }

    //Auxiliary vectors to write HDF5 file
    double *pointVector;
    int *connec2;
    double *pointScalar;
    int *intCellScalar;
    double *douCellScalar;

    pointVector = new double[3*numNodesCoarse]();
    connec2 = new int[6*numElemCoarse]();
    pointScalar = new double[numNodesCoarse]();
    intCellScalar = new int[numElemCoarse]();
    douCellScalar = new double[numElemCoarse]();

    for (int i = 0; i < numNodesCoarse; i++) {
        double *x = nodesCoarse_[i] -> getCoordinates();
        pointVector[3*i  ] = x[0];
        pointVector[3*i+1] = x[1];
        pointVector[3*i+2] = 0.0;
    }
    for (int i = 0; i < numElemCoarse; i++){
        int* con = elementsCoarse_[i] -> getConnectivity();
        connec2[6*i  ] = con[0]; connec2[6*i+1] = con[1];
        connec2[6*i+2] = con[2]; connec2[6*i+3] = con[3];
        connec2[6*i+4] = con[4]; connec2[6*i+5] = con[5];
    }
   
    hid_t file, file2; 
    hid_t dataset; 
    hid_t dataspace;

    herr_t status;
    hsize_t xdim = numNodesCoarse;
    hsize_t pointVectorDims[2] = { numNodesCoarse, 3 };
    hsize_t pointScalarDims[2] = { numNodesCoarse, 1 };
    hsize_t pointCellScalarDims[2] = { numElemCoarse, 1 };
    hsize_t connec2Dims[2] = { numElemCoarse, 6 };

    //Create HDF5 file
    file = H5Fcreate(s1.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (step == 0) file2 = H5Fcreate(s2.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    //Write xdmf file
    output_v << "<?xml version=\"1.0\"?>" << std::endl
             << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>" << std::endl
             << "<Xdmf Version=\"2.0\" xmlns:xi=\"http://www.w3.org/2001/XInclude\" >" << std::endl
             << "<Domain>" << std::endl
             << "  <Grid>"  << std::endl
             << "    <Topology TopologyType=\"Tri_6\" NumberOfElements=\"" << numElemCoarse << "\" >" << std::endl
             << "      <DataItem Format=    \"HDF\" NumberType=\"int\" Dimensions=\"" << numElemCoarse << " 6\" >" << std::endl;
    
    //Connectivity
    output_v << "        " << s2 << ":/connec" << std::endl;
    //Start connectivity
    if (step == 0){
        dataspace = H5Screate_simple(2, connec2Dims, NULL);
        dataset = H5Dcreate2(file2, "/connec", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &connec2[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
    }
    //End connectivity

    output_v << "      </DataItem>" << std::endl
             << "    </Topology>" << std::endl
             << "    <Geometry GeometryType=\"XYZ\">" << std::endl
             << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< numNodesCoarse <<" 3\">" << std::endl;

    //Coordinates   
    output_v << "        " << s2 << ":/coords" << std::endl;
    //Start coordinates
    if (step == 0){
        dataspace = H5Screate_simple(2, pointVectorDims, NULL);
        dataset = H5Dcreate2(file2, "/coords", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointVector[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
    }
    //End coordinates
    output_v << "      </DataItem>" << std::endl
             << "    </Geometry>" << std::endl;

    //LISTS
    //Velocity
    if (coarseModel.printVelocity){
        output_v << "    <Attribute Name=\"Velocity\" Center=\"Node\" AttributeType=\"Vector\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< numNodesCoarse <<" 3\">" << std::endl;  
        output_v << "        " << s1 << ":/velocity" << std::endl;     
        //Start Velocity
        for (int i = 0; i < numNodesCoarse; i++) {
            pointVector[3*i  ] = nodesCoarse_[i] -> getVelocity(0);
            pointVector[3*i+1] = nodesCoarse_[i] -> getVelocity(1);
            pointVector[3*i+2] = 0.0;
        };
        dataspace = H5Screate_simple(2, pointVectorDims, NULL);
        dataset = H5Dcreate2(file, "/velocity", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointVector[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Velocity
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    //Acceleration
    if (coarseModel.printAcceleration){
        output_v << "    <Attribute Name=\"Acceleration\" Center=\"Node\" AttributeType=\"Vector\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< numNodesCoarse <<" 3\">" << std::endl;  
        output_v << "        " << s1 << ":/acceleration" << std::endl;     
        //Start Acceleration
        for (int i = 0; i < numNodesCoarse; i++) {
            pointVector[3*i  ] = nodesCoarse_[i] -> getAcceleration(0);
            pointVector[3*i+1] = nodesCoarse_[i] -> getAcceleration(1);
            pointVector[3*i+2] = 0.0;
        };
        dataspace = H5Screate_simple(2, pointVectorDims, NULL);
        dataset = H5Dcreate2(file, "/acceleration", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointVector[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Acceleration
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    //Real Velocity
    if (coarseModel.printRealVelocity){
        output_v << "    <Attribute Name=\"Real Velocity\" Center=\"Node\" AttributeType=\"Vector\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< numNodesCoarse <<" 3\">" << std::endl;  
        output_v << "        " << s1 << ":/realVelocity" << std::endl;     
        //Start Real Velocity
        for (int i = 0; i < numNodesCoarse; i++) {
            pointVector[3*i  ] = nodesCoarse_[i] -> getVelocityArlequin(0);
            pointVector[3*i+1] = nodesCoarse_[i] -> getVelocityArlequin(1);
            pointVector[3*i+2] = 0.0;
        };
        dataspace = H5Screate_simple(2, pointVectorDims, NULL);
        dataset = H5Dcreate2(file, "/realVelocity", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointVector[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Real Velocity
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    //Lagrange Multipliers
    if (coarseModel.printLagrangeMultipliers){
        output_v << "    <Attribute Name=\"Lagrange Multipliers\" Center=\"Node\" AttributeType=\"Vector\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< numNodesCoarse <<" 3\">" << std::endl;  
        output_v << "        " << s1 << ":/lagrangeMultiplers" << std::endl;     
        //Start Lagrange Multipliers
        for (int i = 0; i < numNodesCoarse; i++) {
            pointVector[3*i  ] = nodesCoarse_[i] -> getLagrangeMultiplier(0);
            pointVector[3*i+1] = nodesCoarse_[i] -> getLagrangeMultiplier(1);
            pointVector[3*i+2] = 0.0;
        };
        dataspace = H5Screate_simple(2, pointVectorDims, NULL);
        dataset = H5Dcreate2(file, "/lagrangeMultiplers", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointVector[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Real Velocity
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    //Dist Function
    if (coarseModel.printDistFunction){
        output_v << "    <Attribute Name=\"Dist Function\" Center=\"Node\" AttributeType=\"Scalar\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\"1 "<< numNodesCoarse <<"\">" << std::endl;  
        output_v << "        " << s1 << ":/distfunction" << std::endl;     
        //Start Dist Function
        for (int i = 0; i < numNodesCoarse; i++) pointScalar[i] = nodesCoarse_[i] -> getDistFunction();
        dataspace = H5Screate_simple(2, pointScalarDims, NULL);
        dataset = H5Dcreate2(file, "/distfunction", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointScalar[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Dist Function
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    //Weight Function
    if (coarseModel.printEnergyWeightFunction){
        output_v << "    <Attribute Name=\"Weight Function\" Center=\"Node\" AttributeType=\"Scalar\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\"1 "<< numNodesCoarse <<"\">" << std::endl;  
        output_v << "        " << s1 << ":/weightFunction" << std::endl;     
        //Start Weigth Function
        for (int i = 0; i < numNodesCoarse; i++) pointScalar[i] = nodesCoarse_[i] -> getWeightFunction();
        dataspace = H5Screate_simple(2, pointScalarDims, NULL);
        dataset = H5Dcreate2(file, "/weightFunction", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointScalar[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Weight Function
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    //Pressure
    if (coarseModel.printPressure){
        output_v << "    <Attribute Name=\"Pressure\" Center=\"Node\" AttributeType=\"Scalar\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< numNodesCoarse <<" 1\">" << std::endl;  
        output_v << "        " << s1 << ":/pressure" << std::endl;     
        //Start Pressure
        for (int i = 0; i < numNodesCoarse; i++) pointScalar[i] = nodesCoarse_[i] -> getPressure();
        dataspace = H5Screate_simple(2, pointScalarDims, NULL);
        dataset = H5Dcreate2(file, "/pressure", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointScalar[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Pressure
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    //Real Pressure
    if (coarseModel.printRealPressure){
        output_v << "    <Attribute Name=\"Real Pressure\" Center=\"Node\" AttributeType=\"Scalar\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< numNodesCoarse <<" 1\">" << std::endl;  
        output_v << "        " << s1 << ":/realPressure" << std::endl;     
        //Start Real Pressure
        for (int i = 0; i < numNodesCoarse; i++) pointScalar[i] = nodesCoarse_[i] -> getPressureArlequin();
        dataspace = H5Screate_simple(2, pointScalarDims, NULL);
        dataset = H5Dcreate2(file, "/realPressure", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointScalar[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Real Pressure
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    if (coarseModel.printVorticity){
        output_v << "    <Attribute Name=\"Vorticity\" Center=\"Node\" AttributeType=\"Scalar\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< numNodesCoarse <<" 1\">" << std::endl;  
        output_v << "        " << s1 << ":/vorticity" << std::endl;     
        //Start Real Pressure
        for (int i = 0; i < numNodesCoarse; i++) pointScalar[i] = nodesCoarse_[i] -> getVorticity();
        dataspace = H5Screate_simple(2, pointScalarDims, NULL);
        dataset = H5Dcreate2(file, "/vorticity", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointScalar[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Real Pressure
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    };

    //Process
    if (coarseModel.printProcess){
        output_v << "    <Attribute Name=\"Process\" Center=\"Cell\" AttributeType=\"Scalar\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"int\" Dimensions=\""<< numElemCoarse <<" 1\">" << std::endl;  
        output_v << "        " << s1 << ":/process" << std::endl;
        //Start Process
        dataspace = H5Screate_simple(2, pointCellScalarDims, NULL);
        dataset = H5Dcreate2(file, "/process", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &domDecompCoarse.first[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Process
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    //WeightFunction2
    if (coarseModel.printEnergyWeightFunction){
        output_v << "    <Attribute Name=\"Weight Function\" Center=\"Cell\" AttributeType=\"Scalar\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< numElemCoarse <<" 1\">" << std::endl;  
        output_v << "        " << s1 << ":/weightFunction2" << std::endl;     
        //Start Weight Function
        for (int i = 0; i < numElemCoarse; i++) douCellScalar[i] = elementsCoarse_[i] -> getIntegPointWeightFunction(0);
        dataspace = H5Screate_simple(2, pointCellScalarDims, NULL);
        dataset = H5Dcreate2(file, "/weightFunction2", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &douCellScalar[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Weight Function
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    //Lines
    if (coarseModel.printLines){
        output_v << "    <Attribute Name=\"Lines\" Center=\"Cell\" AttributeType=\"Scalar\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"int\" Dimensions=\""<< numElemCoarse <<" 1\">" << std::endl;  
        output_v << "        " << s1 << ":/lines" << std::endl;     
        //Start Lines
        for (int i = 0; i < numElemCoarse; i++){
            int res = 0;
            for (int j = 0; j < numBoundElemCoarse; j++)
               if (boundaryCoarse_[j] -> getElement() == i) res = boundaryCoarse_[j] -> getBoundaryGroup();
            intCellScalar[i] = res;
        };
        dataspace = H5Screate_simple(2, pointCellScalarDims, NULL);
        dataset = H5Dcreate2(file, "/lines", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &intCellScalar[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Lines
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    // Jacobian
    if (coarseModel.printJacobian){
        output_v << "    <Attribute Name=\"Jacobian\" Center=\"Cell\" AttributeType=\"Scalar\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< numElemCoarse <<" 1\">" << std::endl;  
        output_v << "        " << s1 << ":/jacobian" << std::endl;     
        //Start Weight Function
        for (int i = 0; i < numElemCoarse; i++) douCellScalar[i] = elementsCoarse_[i] -> getJacobian();
        dataspace = H5Screate_simple(2, pointCellScalarDims, NULL);
        dataset = H5Dcreate2(file, "/jacobian", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &douCellScalar[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Weight Function
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    // END FILE
    output_v << "  </Grid>" << std::endl
             << "</Domain>" << std::endl
             << "</Xdmf>" << std::endl;

    //Delete auxiliary vectors
    delete [] pointVector;
    delete [] connec2;
    delete [] pointScalar;
    delete [] intCellScalar;
    delete [] douCellScalar;

    //End HDF5 file
    status = H5Fclose(file);
    if (step == 0) status = H5Fclose(file2);

    return;
};

//------------------------------------------------------------------------------
//------------------------------PRINT FINE RESULTS------------------------------
//------------------------------------------------------------------------------
template<>
void Arlequin<2,2>::printResultsFine(int step) {

    //PRINT FINE MODEL RESULTS
    std::string result;
    std::ostringstream convert;

    convert << step+100000;
    result = convert.str();

    std::string s = "saidaVelFine"+result+".xdmf";
    std::fstream output_v(s.c_str(), std::ios_base::out);

    std::string s1 = "resultFine"+result+".h5";
    std::fstream filename(s1.c_str(), std::ios_base::out);

    //Auxiliary vectors to write HDF5 file
    double *pointVector;
    int *connec2;
    double *pointScalar;
    int *intCellScalar;
    double *douCellScalar;

    pointVector = new double[3*numNodesFine]();
    connec2 = new int[6*numElemFine]();
    pointScalar = new double[numNodesFine]();
    intCellScalar = new int[numElemFine]();
    douCellScalar = new double[numElemFine]();

    for (int i = 0; i < numNodesFine; i++) {
        double *x = nodesFine_[i] -> getCoordinates();
        pointVector[3*i  ] = x[0];
        pointVector[3*i+1] = x[1];
        pointVector[3*i+2] = 0.0;
    }
    for (int i = 0; i < numElemFine; i++){
        int* con = elementsFine_[i] -> getConnectivity();
        connec2[6*i  ] = con[0]; connec2[6*i+1] = con[1];
        connec2[6*i+2] = con[2]; connec2[6*i+3] = con[3];
        connec2[6*i+4] = con[4]; connec2[6*i+5] = con[5];
    }
   
    hid_t file; 
    hid_t dataset; 
    hid_t dataspace;

    herr_t status;
    hsize_t xdim = numNodesFine;
    hsize_t pointVectorDims[2] = { numNodesFine, 3 };
    hsize_t pointScalarDims[2] = { numNodesFine, 1 };
    hsize_t pointCellScalarDims[2] = { numElemFine, 1 };
    hsize_t connec2Dims[2] = { numElemFine, 6 };

    //Create HDF5 file
    file = H5Fcreate(s1.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    //Write xdmf file
    output_v << "<?xml version=\"1.0\"?>" << std::endl
             << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>" << std::endl
             << "<Xdmf Version=\"2.0\" xmlns:xi=\"http://www.w3.org/2001/XInclude\" >" << std::endl
             << "<Domain>" << std::endl
             << "  <Grid>"  << std::endl
             << "    <Topology TopologyType=\"Tri_6\" NumberOfElements=\"" << numElemFine << "\" >" << std::endl
             << "      <DataItem Format=    \"HDF\" NumberType=\"int\" Dimensions=\"" << numElemFine << " 6\" >" << std::endl;
    
    //Connectivity
    output_v << "        " << s1 << ":/connec" << std::endl;
    //Start connectivity
    dataspace = H5Screate_simple(2, connec2Dims, NULL);
    dataset = H5Dcreate2(file, "/connec", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &connec2[0]);
    status = H5Dclose(dataset);
    status = H5Sclose(dataspace);
    //End connectivity

    output_v << "      </DataItem>" << std::endl
             << "    </Topology>" << std::endl
             << "    <Geometry GeometryType=\"XYZ\">" << std::endl
             << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< numNodesFine <<" 3\">" << std::endl;

    //Coordinates   
    output_v << "        " << s1 << ":/coords" << std::endl;
    //Start coordinates
    dataspace = H5Screate_simple(2, pointVectorDims, NULL);
    dataset = H5Dcreate2(file, "/coords", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointVector[0]);
    status = H5Dclose(dataset);
    status = H5Sclose(dataspace);
    //End coordinates
    output_v << "      </DataItem>" << std::endl
             << "    </Geometry>" << std::endl;

    //LISTS
    //Velocity
    if (fineModel.printVelocity){
        output_v << "    <Attribute Name=\"Velocity\" Center=\"Node\" AttributeType=\"Vector\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< numNodesFine <<" 3\">" << std::endl;  
        output_v << "        " << s1 << ":/velocity" << std::endl;     
        //Start Velocity
        for (int i = 0; i < numNodesFine; i++) {
            pointVector[3*i  ] = nodesFine_[i] -> getVelocity(0);
            pointVector[3*i+1] = nodesFine_[i] -> getVelocity(1);
            pointVector[3*i+2] = 0.0;
        };
        dataspace = H5Screate_simple(2, pointVectorDims, NULL);
        dataset = H5Dcreate2(file, "/velocity", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointVector[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Velocity
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    //Acceleration
    if (fineModel.printAcceleration){
        output_v << "    <Attribute Name=\"Acceleration\" Center=\"Node\" AttributeType=\"Vector\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< numNodesFine <<" 3\">" << std::endl;  
        output_v << "        " << s1 << ":/acceleration" << std::endl;     
        //Start Acceleration
        for (int i = 0; i < numNodesFine; i++) {
            pointVector[3*i  ] = nodesFine_[i] -> getAcceleration(0);
            pointVector[3*i+1] = nodesFine_[i] -> getAcceleration(1);
            pointVector[3*i+2] = 0.0;
        };
        dataspace = H5Screate_simple(2, pointVectorDims, NULL);
        dataset = H5Dcreate2(file, "/acceleration", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointVector[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Acceleration
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    //Real Velocity
    if (fineModel.printRealVelocity){
        output_v << "    <Attribute Name=\"Real Velocity\" Center=\"Node\" AttributeType=\"Vector\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< numNodesFine <<" 3\">" << std::endl;  
        output_v << "        " << s1 << ":/realVelocity" << std::endl;     
        //Start Real Velocity
        for (int i = 0; i < numNodesFine; i++) {
            pointVector[3*i  ] = nodesFine_[i] -> getVelocityArlequin(0);
            pointVector[3*i+1] = nodesFine_[i] -> getVelocityArlequin(1);
            pointVector[3*i+2] = 0.0;
        };
        dataspace = H5Screate_simple(2, pointVectorDims, NULL);
        dataset = H5Dcreate2(file, "/realVelocity", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointVector[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Real Velocity
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    //Lagrange Multipliers
    if (fineModel.printLagrangeMultipliers){
        output_v << "    <Attribute Name=\"Lagrange Multipliers\" Center=\"Node\" AttributeType=\"Vector\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< numNodesFine <<" 3\">" << std::endl;  
        output_v << "        " << s1 << ":/lagrangeMultiplers" << std::endl;     
        //Start Lagrange Multipliers
        for (int i = 0; i < numNodesFine; i++) {
            pointVector[3*i  ] = nodesFine_[i] -> getLagrangeMultiplier(0);
            pointVector[3*i+1] = nodesFine_[i] -> getLagrangeMultiplier(1);
            pointVector[3*i+2] = 0.0;
        };
        dataspace = H5Screate_simple(2, pointVectorDims, NULL);
        dataset = H5Dcreate2(file, "/lagrangeMultiplers", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointVector[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Lagrange Multipliers
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    //Normal vector
    if (fineModel.printInnerNormal){
        output_v << "    <Attribute Name=\"Normal Vector\" Center=\"Node\" AttributeType=\"Vector\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< numNodesFine <<" 3\">" << std::endl;  
        output_v << "        " << s1 << ":/normalVector" << std::endl;     
        //Start Normal vector
        for (int i = 0; i < numNodesFine; i++) {
            double *n = nodesFine_[i] -> getInnerNormal();
            pointVector[3*i  ] = n[0];
            pointVector[3*i+1] = n[1];
            pointVector[3*i+2] = 0.0;
        };
        dataspace = H5Screate_simple(2, pointVectorDims, NULL);
        dataset = H5Dcreate2(file, "/normalVector", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointVector[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Normal Vector
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    //Mesh Velocity
    if (fineModel.printMeshVelocity){
        output_v << "    <Attribute Name=\"Mesh Velocity\" Center=\"Node\" AttributeType=\"Vector\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< numNodesFine <<" 3\">" << std::endl;  
        output_v << "        " << s1 << ":/meshVelocity" << std::endl;     
        //Start mesh Velocity
        for (int i = 0; i < numNodesFine; i++) {
            pointVector[3*i  ] = nodesFine_[i] -> getMeshVelocity(0);
            pointVector[3*i+1] = nodesFine_[i] -> getMeshVelocity(1);
            pointVector[3*i+2] = 0.0;
        };
        dataspace = H5Screate_simple(2, pointVectorDims, NULL);
        dataset = H5Dcreate2(file, "/meshVelocity", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointVector[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Mesh Velocity
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    };

    //Mesh Displacement
    if (fineModel.printMeshDisplacement){
        output_v << "    <Attribute Name=\"Mesh Displacement\" Center=\"Node\" AttributeType=\"Vector\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< numNodesFine <<" 3\">" << std::endl;  
        output_v << "        " << s1 << ":/meshDisplacement" << std::endl;     
        //Start Mesh displacement
        for (int i = 0; i < numNodesFine; i++) {
            double* x = nodesFine_[i] -> getCoordinates();
            double *xi = nodesFine_[i] -> getInitialCoordinates();
            pointVector[3*i  ] = x[0] - xi[0];
            pointVector[3*i+1] = x[1] - xi[1];
            pointVector[3*i+2] = 0.0;
        };
        dataspace = H5Screate_simple(2, pointVectorDims, NULL);
        dataset = H5Dcreate2(file, "/meshDisplacement", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointVector[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Mesh displacement
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    };

    //Element correspondence
    if (fineModel.printElementCorrespondence){
        output_v << "    <Attribute Name=\"Element\" Center=\"Node\" AttributeType=\"Scalar\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\"1 "<< numNodesFine <<"\">" << std::endl;  
        output_v << "        " << s1 << ":/elementCorresp" << std::endl;     
        //Start Element correspondence
        for (int i = 0; i < numNodesFine; i++) pointScalar[i] = nodesFine_[i] -> getNodalElemCorrespondence();
        dataspace = H5Screate_simple(2, pointScalarDims, NULL);
        dataset = H5Dcreate2(file, "/elementCorresp", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointScalar[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Element correspondence
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    //Dist Function
    if (fineModel.printDistFunction){
        output_v << "    <Attribute Name=\"Dist Function\" Center=\"Node\" AttributeType=\"Scalar\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\"1 "<< numNodesFine <<"\">" << std::endl;  
        output_v << "        " << s1 << ":/distfunction" << std::endl;     
        //Start Dist Function
        for (int i = 0; i < numNodesFine; i++) pointScalar[i] = nodesFine_[i] -> getDistFunction();
        dataspace = H5Screate_simple(2, pointScalarDims, NULL);
        dataset = H5Dcreate2(file, "/distfunction", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointScalar[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Dist Function
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    //Weight Function
    if (fineModel.printEnergyWeightFunction){
        output_v << "    <Attribute Name=\"Weight Function\" Center=\"Node\" AttributeType=\"Scalar\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\"1 "<< numNodesFine <<"\">" << std::endl;  
        output_v << "        " << s1 << ":/weightFunction" << std::endl;     
        //Start Weigth Function
        for (int i = 0; i < numNodesFine; i++) pointScalar[i] = nodesFine_[i] -> getWeightFunction();
        dataspace = H5Screate_simple(2, pointScalarDims, NULL);
        dataset = H5Dcreate2(file, "/weightFunction", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointScalar[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Weight Function
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    //Pressure
    if (fineModel.printPressure){
        output_v << "    <Attribute Name=\"Pressure\" Center=\"Node\" AttributeType=\"Scalar\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< numNodesFine <<" 1\">" << std::endl;  
        output_v << "        " << s1 << ":/pressure" << std::endl;     
        //Start Pressure
        for (int i = 0; i < numNodesFine; i++) pointScalar[i] = nodesFine_[i] -> getPressure();
        dataspace = H5Screate_simple(2, pointScalarDims, NULL);
        dataset = H5Dcreate2(file, "/pressure", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointScalar[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Pressure
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    //Real Pressure
    if (fineModel.printRealPressure){
        output_v << "    <Attribute Name=\"Real Pressure\" Center=\"Node\" AttributeType=\"Scalar\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< numNodesFine <<" 1\">" << std::endl;  
        output_v << "        " << s1 << ":/realPressure" << std::endl;     
        //Start Real Pressure
        for (int i = 0; i < numNodesFine; i++) pointScalar[i] = nodesFine_[i] -> getPressureArlequin();
        dataspace = H5Screate_simple(2, pointScalarDims, NULL);
        dataset = H5Dcreate2(file, "/realPressure", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointScalar[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Real Pressure
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    if (fineModel.printVorticity){
        output_v << "    <Attribute Name=\"Vorticity\" Center=\"Node\" AttributeType=\"Scalar\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< numNodesFine <<" 1\">" << std::endl;  
        output_v << "        " << s1 << ":/vorticity" << std::endl;     
        //Start Real Pressure
        for (int i = 0; i < numNodesFine; i++) pointScalar[i] = nodesFine_[i] -> getVorticity();
        dataspace = H5Screate_simple(2, pointScalarDims, NULL);
        dataset = H5Dcreate2(file, "/vorticity", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointScalar[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Real Pressure
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    };

    //Process
    if (fineModel.printProcess){
        output_v << "    <Attribute Name=\"Process\" Center=\"Cell\" AttributeType=\"Scalar\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"int\" Dimensions=\""<< numElemFine <<" 1\">" << std::endl;  
        output_v << "        " << s1 << ":/process" << std::endl;
        //Start Process
        dataspace = H5Screate_simple(2, pointCellScalarDims, NULL);
        dataset = H5Dcreate2(file, "/process", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &domDecompFine.first[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Process
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    //WeightFunction2
    if (fineModel.printEnergyWeightFunction){
        output_v << "    <Attribute Name=\"Weight Function\" Center=\"Cell\" AttributeType=\"Scalar\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< numElemFine <<" 1\">" << std::endl;  
        output_v << "        " << s1 << ":/weightFunction2" << std::endl;     
        //Start Weight Function
        for (int i = 0; i < numElemFine; i++) douCellScalar[i] = elementsFine_[i] -> getIntegPointWeightFunction(0);
        dataspace = H5Screate_simple(2, pointCellScalarDims, NULL);
        dataset = H5Dcreate2(file, "/weightFunction2", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &douCellScalar[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Weight Function
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    //Lines
    if (fineModel.printLines){
        output_v << "    <Attribute Name=\"Lines\" Center=\"Cell\" AttributeType=\"Scalar\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"int\" Dimensions=\""<< numElemFine <<" 1\">" << std::endl;  
        output_v << "        " << s1 << ":/lines" << std::endl;     
        //Start Lines
        for (int i = 0; i < numElemFine; i++){
            int res = 0;
            for (int j = 0; j < numBoundElemFine; j++)
               if (boundaryFine_[j] -> getElement() == i) res = boundaryFine_[j] -> getBoundaryGroup();
            intCellScalar[i] = res;
        };
        dataspace = H5Screate_simple(2, pointCellScalarDims, NULL);
        dataset = H5Dcreate2(file, "/lines", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &intCellScalar[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Lines
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    // Jacobian
    if (fineModel.printJacobian){
        output_v << "    <Attribute Name=\"Jacobian\" Center=\"Cell\" AttributeType=\"Scalar\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< numElemFine <<" 1\">" << std::endl;  
        output_v << "        " << s1 << ":/jacobian" << std::endl;     
        //Start Weight Function
        for (int i = 0; i < numElemFine; i++) douCellScalar[i] = elementsFine_[i] -> getJacobian();
        dataspace = H5Screate_simple(2, pointCellScalarDims, NULL);
        dataset = H5Dcreate2(file, "/jacobian", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &douCellScalar[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Weight Function
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    //Glue Zone
    if (fineModel.printGlueZone){
        output_v << "    <Attribute Name=\"Glue Zone\" Center=\"Cell\" AttributeType=\"Scalar\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"int\" Dimensions=\""<< numElemFine <<" 1\">" << std::endl;  
        output_v << "        " << s1 << ":/glueZone" << std::endl;     
        //Start Lines
        for (int i = 0; i < numElemFine; i++) intCellScalar[i] = elementsFine_[i] -> getGlueZoneInt();
        dataspace = H5Screate_simple(2, pointCellScalarDims, NULL);
        dataset = H5Dcreate2(file, "/glueZone", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &intCellScalar[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Lines
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    // END FILE
    output_v << "  </Grid>" << std::endl
             << "</Domain>" << std::endl
             << "</Xdmf>" << std::endl;

    //Delete auxiliary vectors
    delete [] pointVector;
    delete [] connec2;
    delete [] pointScalar;
    delete [] intCellScalar;
    delete [] douCellScalar;

    //End HDF5 file
    status = H5Fclose(file);

    return;
};


//------------------------------------------------------------------------------
//------------SETS COARSE/FINE MESHES AND GETS ITS BASIC INFORMATIONS-----------
//------------------------------------------------------------------------------
template<>
void Arlequin<2,2>::setFluidModels(FluidMesh& coarse, FluidMesh& fine){

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

    // // Rotating fine mesh
    // for (int i = 0; i < numNodesFine; ++i){
    //     typename Nodes::VecLocD x, xn;
    //     x = nodesFine_[i] -> getCoordinates();       
    
    //     double a = -90 * pi / 180;

    //     xn(0) = 0.0 + (x(0)-0.0) * cos(a) - (x(1)-0.0) * sin(a);
    //     xn(1) = 0.0 + (x(0)-0.0) * sin(a) + (x(1)-0.0) * cos(a);

    //     nodesFine_[i] -> setCoordinates(xn);
    //     nodesFine_[i] -> setPreviousCoordinates(0,xn(0));
    //     nodesFine_[i] -> setPreviousCoordinates(1,xn(1));
    // }

    parametersFine = &fineModel.fluidParameters;
    parametersCoarse = &coarseModel.fluidParameters;

    setSignaledDistance();

    //Construct the glue zone based on some defined criterion
    setCouplingZone();

    //Computes the Weight function for all the finite elements
    setWeightFunction(16.); 

    //Update domain decomposition - Start
    int size;
    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    Vec  b;
    PetscErrorCode    ierr;
    int start[size], end[size];
    int numDOF = 3*numNodesCoarse + 3*numNodesFine + 2*numNodesGlueZoneFine;

    ierr = VecCreate(PETSC_COMM_WORLD,&b);
    ierr = VecSetSizes(b,PETSC_DECIDE,numDOF);
    VecSetFromOptions(b);
    ierr = VecGetOwnershipRange(b,&start[rank],&end[rank]);

    VecDestroy(&b);
    
    for (int i = 0; i < size; ++i){
        MPI_Bcast(&start[i],1,MPI_INT,i,PETSC_COMM_WORLD);
        MPI_Bcast(&end[i],1,MPI_INT,i,PETSC_COMM_WORLD);
    }

    for (int i = 0; i < numElemCoarse; i++){
        int* connec = elementsCoarse_[i] -> getConnectivity();
        for (int j = 0; j < size; j++){
            if ((3*connec[3] >= start[j]) && (3*connec[3] <= end[j])) {
                domDecompCoarse.first[i] = j;
                break;
            }
        }
    }    
    for (int i = 0; i < numElemFine; i++){
        int* connec = elementsFine_[i] -> getConnectivity();
        for (int j = 0; j < size; j++){
            if ((3*(numNodesCoarse + connec[3]) >= start[j]) && (3*(numNodesCoarse + connec[3]) <= end[j])) {
                domDecompFine.first[i] = j;
                break;
            }
        }
    }    

    //Update domain decomposition - End



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
void Arlequin<2,2>::dragAndLiftCoefficients(std::ofstream& dragLift){



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
        double aux_Mom = 0.;
        double aux_Per = 0.;



        for (int i=0; i<fineModel.numberOfLines; i++){
            if (boundaryFine_[jel] -> getBoundaryGroup() == fineModel.dragAndLiftBoundary[i]){
                
                int iel = boundaryFine_[jel] -> getElement();
                elementsFine_[iel] -> computeDragAndLiftForces(pDForce, pLForce, fDForce, fLForce, dForce, lForce, aux_Mom, aux_Per);

                pMom += aux_Mom;
                per += aux_Per;
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

    return;
}

//------------------------------------------------------------------------------
//----------------COMPUTE ARLEQUIN COUPLED NAVIER-STOKES PROBLEM----------------
//------------------------------------------------------------------------------
template<>
void Arlequin<2,2>::setMatVecValuesCoarseModel(double **matrix, double *rhs, int* connec){

    //Disperse local contributions into the global matrix
    for (int i = 0; i < 6; i++){
        for (int j = 0; j < 6; j++){
            //Matrix K and C            
            int dof_i = 3 * connec[i];
            int dof_j = 3 * connec[j];
            MatSetValues(A,1,&dof_i,1,&dof_j,&matrix[2*i  ][2*j  ],ADD_VALUES);
            
            dof_i = 3 * connec[i] + 1;
            dof_j = 3 * connec[j];
            MatSetValues(A,1,&dof_i,1,&dof_j,&matrix[2*i+1][2*j  ],ADD_VALUES);
            
            dof_i = 3 * connec[i];
            dof_j = 3 * connec[j] + 1;
            MatSetValues(A,1,&dof_i,1,&dof_j,&matrix[2*i  ][2*j+1],ADD_VALUES);
            
            dof_i = 3 * connec[i] + 1;
            dof_j = 3 * connec[j] + 1;
            MatSetValues(A,1,&dof_i,1,&dof_j,&matrix[2*i+1][2*j+1],ADD_VALUES);
        
            //Matrix G and Gt
            dof_i = 3 * connec[i];
            dof_j = 3 * connec[j] + 2;
            MatSetValues(A,1,&dof_i,1,&dof_j,&matrix[2*i  ][12+j],ADD_VALUES);
            MatSetValues(A,1,&dof_j,1,&dof_i,&matrix[12+j][2*i  ],ADD_VALUES);
            
            dof_i = 3 * connec[i] + 1;
            dof_j = 3 * connec[j] + 2;
            MatSetValues(A,1,&dof_i,1,&dof_j,&matrix[2*i+1][12+j],ADD_VALUES);
            MatSetValues(A,1,&dof_j,1,&dof_i,&matrix[12+j][2*i+1],ADD_VALUES);

            // Matrix Q
            dof_i = 3 * connec[i] + 2;
            dof_j = 3 * connec[j] + 2;
            MatSetValues(A,1,&dof_j,1,&dof_i,&matrix[12+i][12+j],ADD_VALUES);
        };                  
        //Rhs vector
        int dof_i = 3 * connec[i];
        VecSetValues(b,1,&dof_i,&rhs[2*i  ],ADD_VALUES);
        
        dof_i = 3 * connec[i] + 1;
        VecSetValues(b,1,&dof_i,&rhs[2*i+1],ADD_VALUES);

        dof_i = 3 * connec[i] + 2;
        VecSetValues(b,1,&dof_i,&rhs[12+i],ADD_VALUES);
    };

    return;
}

//------------------------------------------------------------------------------
//----------------COMPUTE ARLEQUIN COUPLED NAVIER-STOKES PROBLEM----------------
//------------------------------------------------------------------------------
template<>
void Arlequin<2,2>::setMatVecValuesFineModel(double **matrix, double *rhs, int* connec){

    //Disperse local contributions into the global matrix
    for (int i = 0; i < 6; i++){
        for (int j = 0; j < 6; j++){
            //Matrix K and C
            int dof_i = 3*numNodesCoarse + 3 * connec[i];
            int dof_j = 3*numNodesCoarse + 3 * connec[j];
            MatSetValues(A,1,&dof_i,1,&dof_j,&matrix[2*i  ][2*j  ],ADD_VALUES);
            
            dof_i = 3*numNodesCoarse + 3 * connec[i] + 1;
            dof_j = 3*numNodesCoarse + 3 * connec[j];
            MatSetValues(A,1,&dof_i,1,&dof_j,&matrix[2*i+1][2*j  ],ADD_VALUES);
            
            dof_i = 3*numNodesCoarse + 3 * connec[i];
            dof_j = 3*numNodesCoarse + 3 * connec[j] + 1;
            MatSetValues(A,1,&dof_i,1,&dof_j,&matrix[2*i  ][2*j+1],ADD_VALUES);

            dof_i = 3*numNodesCoarse + 3 * connec[i] + 1;
            dof_j = 3*numNodesCoarse + 3 * connec[j] + 1;
            MatSetValues(A,1,&dof_i,1,&dof_j,&matrix[2*i+1][2*j+1],ADD_VALUES);
            
            //Matrix G and Gt
            dof_i = 3 * numNodesCoarse + 3 * connec[i];
            dof_j = 3 * numNodesCoarse + 3 * connec[j] + 2;
            MatSetValues(A,1,&dof_i,1,&dof_j,&matrix[2*i  ][12+j],ADD_VALUES);
            MatSetValues(A,1,&dof_j,1,&dof_i,&matrix[12+j][2*i  ],ADD_VALUES);
            
            dof_i = 3 * numNodesCoarse + 3 * connec[i] + 1;
            dof_j = 3 * numNodesCoarse + 3 * connec[j] + 2;
            MatSetValues(A,1,&dof_i,1,&dof_j,&matrix[2*i+1][12+j],ADD_VALUES);
            MatSetValues(A,1,&dof_j,1,&dof_i,&matrix[12+j][2*i+1],ADD_VALUES);
            
            //Matrix Q
            dof_i = 3 * numNodesCoarse + 3 * connec[i] + 2;
            dof_j = 3 * numNodesCoarse + 3 * connec[j] + 2;
            MatSetValues(A,1,&dof_j,1,&dof_i,&matrix[12+i][12+j],ADD_VALUES);
        };     
        ///Rhs vector
        int dof_i = 3 * numNodesCoarse + 3 * connec[i];
        VecSetValues(b,1,&dof_i,&rhs[2*i  ],ADD_VALUES);
        
        dof_i = 3 * numNodesCoarse + 3 * connec[i] + 1;
        VecSetValues(b,1,&dof_i,&rhs[2*i+1],ADD_VALUES);

        dof_i = 3 * numNodesCoarse + 3 * connec[i] + 2;
        VecSetValues(b,1,&dof_i,&rhs[12+i],ADD_VALUES);
    }; 

    return;
};

//------------------------------------------------------------------------------
//----------------COMPUTE ARLEQUIN COUPLED NAVIER-STOKES PROBLEM----------------
//------------------------------------------------------------------------------
template<>
void Arlequin<2,2>::setMatVecValuesLagMultFineFine(double **Ajac2, double **localMV_mat, 
                                                 double **ArlequinA1, double **ArlequinA2, 
                                                 double *Rhs2, double *rhsLagMult2,
                                                 double *localMV_vec, double *RhsArlequin2,
                                                 int* connec, int* connecL){

    double &alpha_f = parametersFine -> getAlphaF();
    double &alpha_m = parametersFine -> getAlphaM();
    double &gamma = parametersFine -> getGamma();
    
    double integ = alpha_f * gamma * dTime;
    //Disperse local contributions into the global matrix
    for (int i=0; i<6; i++){
        for (int j=0; j<6; j++){
            //COUPLING OPERATOR
            if (fabs(Ajac2[2*i  ][2*j  ]) >= 1.e-15){
                int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL[i];
                int d_j = 3*numNodesCoarse + 3*connec[j];
                double value = Ajac2[2*i  ][2*j  ]*integ;
                MatSetValues(A,1,&d_i,1,&d_j,&value,ADD_VALUES);
                MatSetValues(A,1,&d_j,1,&d_i,&Ajac2[2*i  ][2*j  ],ADD_VALUES);
            };
            if (fabs(Ajac2[2*i+1][2*j  ]) >= 1.e-15){
                int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL[i] + 1;
                int d_j = 3*numNodesCoarse + 3*connec[j];
                double value = Ajac2[2*i+1][2*j  ]*integ;
                MatSetValues(A,1,&d_i,1,&d_j,&value,ADD_VALUES);
                MatSetValues(A,1,&d_j,1,&d_i,&Ajac2[2*i+1][2*j  ],ADD_VALUES);
            };
            if (fabs(Ajac2[2*i+1][2*j+1]) >= 1.e-15){
                int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL[i] + 1;
                int d_j = 3*numNodesCoarse + 3*connec[j] + 1;
                double value = Ajac2[2*i+1][2*j+1]*integ;
                MatSetValues(A,1,&d_i,1,&d_j,&value,ADD_VALUES);
                MatSetValues(A,1,&d_j,1,&d_i,&Ajac2[2*i+1][2*j+1],ADD_VALUES);
            };
            if (fabs(Ajac2[2*i  ][2*j+1]) >= 1.e-15){
                int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL[i];
                int d_j = 3*numNodesCoarse + 3*connec[j] + 1;
                double value = Ajac2[2*i  ][2*j+1]*integ;
                MatSetValues(A,1,&d_i,1,&d_j,&value,ADD_VALUES);
                MatSetValues(A,1,&d_j,1,&d_i,&Ajac2[2*i  ][2*j+1],ADD_VALUES);
            };

            //SUPG STABILIZATION
            if (fabs(localMV_mat[2*i  ][2*j  ]) >= 1.e-15){
                int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL[i];
                int d_j = 3*numNodesCoarse + 3*connec[j];
                MatSetValues(A,1,&d_j,1,&d_i,&localMV_mat[2*i  ][2*j  ],ADD_VALUES);
            };
            if (fabs(localMV_mat[2*i+1][2*j+1]) >= 1.e-15){
                int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL[i] + 1;
                int d_j = 3*numNodesCoarse + 3*connec[j] + 1;
                MatSetValues(A,1,&d_j,1,&d_i,&localMV_mat[2*i+1][2*j+1],ADD_VALUES);
            };

            //PSPG STABILIZATION
            if (fabs(localMV_mat[2*i  ][12+j]) >= 1.e-15){
                int dof_i = 3*numNodesCoarse + 3*connec[i] + 2;
                int dof_j = 3*numNodesCoarse + 3*numNodesFine + connecL[j];
                MatSetValues(A,1,&dof_i,1,&dof_j,&localMV_mat[2*i  ][12+j],ADD_VALUES);
            };
            if (fabs(localMV_mat[2*i+1][12+j]) >= 1.e-15){
                int dof_i = 3*numNodesCoarse + 3*connec[i] + 2;
                int dof_j = 3*numNodesCoarse + 3*numNodesFine + connecL[j] + 1;
                MatSetValues(A,1,&dof_i,1,&dof_j,&localMV_mat[2*i+1][12+j],ADD_VALUES);
            };

            //ARLEQUIN STABILIZATION
            if (fabs(ArlequinA2[2*i  ][2*j  ]) >= 1.e-15){
                int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL[i];
                int d_j = 3*numNodesCoarse + 3*connec[j];
                MatSetValues(A,1,&d_i,1,&d_j,&ArlequinA2[2*i  ][2*j  ],ADD_VALUES);
            };
            if (fabs(ArlequinA2[2*i+1][2*j+1]) >= 1.e-15){
                int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL[i] + 1;
                int d_j = 3*numNodesCoarse + 3*connec[j] + 1;
                MatSetValues(A,1,&d_i,1,&d_j,&ArlequinA2[2*i+1][2*j+1],ADD_VALUES);
            };
            if (fabs(ArlequinA2[12+i][2*j  ]) >= 1.e-15){
                int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL[i];
                int dof_j = 3*numNodesCoarse + 3*connec[j] + 2;
                MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinA2[12+i][2*j  ],ADD_VALUES);
            };
            if (fabs(ArlequinA2[12+i][2*j+1]) >= 1.e-15){
                int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL[i] + 1;
                int dof_j = 3*numNodesCoarse + 3*connec[j] + 2;
                MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinA2[12+i][2*j+1],ADD_VALUES);
            };
            if (fabs(ArlequinA1[2*i  ][2*j  ]) >= 1.e-15){
                int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL[i];
                int dof_j = 3*numNodesCoarse + 3*numNodesFine + 2*connecL[j];
                MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinA1[2*i  ][2*j  ],ADD_VALUES);
                dof_i++; dof_j++;
                MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinA1[2*i+1][2*j+1],ADD_VALUES);
            };                            
        };

        //RHS VECTOR
        //COUPLING OPERATOR
        int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL[i];
        VecSetValues(b,1,&dof_i,&Rhs2[2*i  ],ADD_VALUES);

        dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL[i] + 1;
        VecSetValues(b,1,&dof_i,&Rhs2[2*i+1],ADD_VALUES);

        dof_i = 3*numNodesCoarse + 3*connec[i];
        VecSetValues(b,1,&dof_i,&rhsLagMult2[2*i  ],ADD_VALUES);

        dof_i = 3*numNodesCoarse + 3*connec[i] + 1;
        VecSetValues(b,1,&dof_i,&rhsLagMult2[2*i+1],ADD_VALUES);

        //SUPG STABILIZATION
        dof_i = 3*numNodesCoarse + 3*connec[i];                       
        VecSetValues(b,1,&dof_i,&localMV_vec[2*i  ],ADD_VALUES);

        dof_i = 3*numNodesCoarse + 3*connec[i] + 1;
        VecSetValues(b,1,&dof_i,&localMV_vec[2*i+1],ADD_VALUES);

        //PSPG STABILIZATION
        dof_i = 3*numNodesCoarse + 3*connec[i] + 2;
        VecSetValues(b,1,&dof_i,&localMV_vec[12+i],ADD_VALUES);

        // ///ARLEQUIN STABILIZATION
        dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL[i];
        VecSetValues(b,1,&dof_i,&RhsArlequin2[2*i  ],ADD_VALUES);
        dof_i++;
        VecSetValues(b,1,&dof_i,&RhsArlequin2[2*i+1],ADD_VALUES);
    };


    return;
};

//------------------------------------------------------------------------------
//----------------COMPUTE ARLEQUIN COUPLED NAVIER-STOKES PROBLEM----------------
//------------------------------------------------------------------------------
template<>
void Arlequin<2,2>::setMatVecValuesLagMultFineCoarse(double **Ajac2, double **localMV_mat, 
                                                   double **ArlequinA1, double **ArlequinA2, 
                                                   double *Rhs2, double *rhsLagMult2,
                                                   double *localMV_vec, double *RhsArlequin2,
                                                   int* connecC, int* connecL){

    double &alpha_f = parametersFine -> getAlphaF();
    double &alpha_m = parametersFine -> getAlphaM();
    double &gamma = parametersFine -> getGamma();
    
    double integ = alpha_f * gamma * dTime;

    //Disperse local contribution into the global matrix
    for (int i=0; i<6; i++){
        for (int j=0; j<6; j++){
            //COUPLING OPERATOR
            if (fabs(Ajac2[2*i  ][2*j  ]) >= 1.e-15){
                int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL[i];
                int dof_j = 3*connecC[j];
                double value = Ajac2[2*i  ][2*j  ] * integ;
                MatSetValues(A,1,&dof_i,1,&dof_j,&value,ADD_VALUES);
                MatSetValues(A,1,&dof_j,1,&dof_i,&Ajac2[2*i  ][2*j  ],ADD_VALUES);
            }
            if (fabs(Ajac2[2*i+1][2*j  ]) >= 1.e-15){
                int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL[i] + 1;
                int dof_j = 3*connecC[j];
                double value = Ajac2[2*i+1][2*j  ]*integ;
                MatSetValues(A,1,&dof_i,1,&dof_j,&value,ADD_VALUES);
                MatSetValues(A,1,&dof_j,1,&dof_i,&Ajac2[2*i+1][2*j  ],ADD_VALUES);
            };
            if (fabs(Ajac2[2*i+1][2*j+1]) >= 1.e-15){
                int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL[i] + 1;
                int dof_j = 3*connecC[j] + 1;
                double value = Ajac2[2*i+1][2*j+1]*integ;
                MatSetValues(A,1,&dof_i,1,&dof_j,&value,ADD_VALUES);
                MatSetValues(A,1,&dof_j,1,&dof_i,&Ajac2[2*i+1][2*j+1],ADD_VALUES);
            };
            if (fabs(Ajac2[2*i  ][2*j+1]) >= 1.e-15){
                int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL[i];
                int dof_j = 3*connecC[j] + 1;
                double value = Ajac2[2*i  ][2*j+1]*integ;
                MatSetValues(A,1,&dof_i,1,&dof_j,&value,ADD_VALUES);
                MatSetValues(A,1,&dof_j,1,&dof_i,&Ajac2[2*i  ][2*j+1],ADD_VALUES);
            };
            // if (fabs(Ajac(12+i,12+j)) >= 1.e-15){
            //     int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL[i];
            //     int dof_j = 3*numNodesCoarse + 3*numNodesFine + 2*connecL[j];
            //     dof_i++; dof_j++;
            //     ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&Ajac(12+i,12+j),ADD_VALUES);
            // };
     

            //SUPG STABILIZATION
            if (fabs(localMV_mat[2*i  ][2*j  ]) >= 1.e-15){
                int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL[i];
                int dof_j = 3*connecC[j];
                MatSetValues(A,1,&dof_j,1,&dof_i,&localMV_mat[2*i  ][2*j  ],ADD_VALUES);
            };
            if (fabs(localMV_mat[2*i+1][2*j+1]) >= 1.e-15){
                int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL[i] + 1;
                int dof_j = 3*connecC[j] + 1;
                MatSetValues(A,1,&dof_j,1,&dof_i,&localMV_mat[2*i+1][2*j+1],ADD_VALUES);
            };

            //PSPG STABILIZATION
            if (fabs(localMV_mat[2*i  ][12+j]) >= 1.e-15){
                int dof_i = 3*connecC[i] + 2;
                int dof_j = 3*numNodesCoarse + 2*numNodesFine + connecL[j];
                MatSetValues(A,1,&dof_i,1,&dof_j,&localMV_mat[2*i  ][12+j],ADD_VALUES);
            };
            if (fabs(localMV_mat[2*i+1][12+j]) >= 1.e-15){
                int dof_i = 3*connecC[i] + 2;
                int dof_j = 3*numNodesCoarse + 2*numNodesFine + connecL[j] + 1;
                MatSetValues(A,1,&dof_i,1,&dof_j,&localMV_mat[2*i+1][12+j],ADD_VALUES);
            };

            //ARLEQUIN STABILIZATION
            if (fabs(ArlequinA2[2*i  ][2*j  ]) >= 1.e-15){
                int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL[i];
                int d_j = 3*connecC[j];
                MatSetValues(A,1,&d_i,1,&d_j,&ArlequinA2[2*i  ][2*j  ],ADD_VALUES);
            };
            if (fabs(ArlequinA2[2*i+1][2*j+1]) >= 1.e-15){
                int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL[i] + 1;
                int d_j = 3*connecC[j] + 1;
                MatSetValues(A,1,&d_i,1,&d_j,&ArlequinA2[2*i+1][2*j+1],ADD_VALUES);
            };

            if (fabs(ArlequinA2[12+i][2*j  ]) >= 1.e-15){
                int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL[i];
                int dof_j = 3*connecC[j] + 2;
                MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinA2[12+i][2*j  ],ADD_VALUES);
            };
            if (fabs(ArlequinA2[12+i][2*j+1]) >= 1.e-15){
                int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL[i] + 1;
                int dof_j = 2*connecC[j] + 2;
                MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinA2[12+i][2*j+1],ADD_VALUES);
            };

            if (fabs(ArlequinA1[2*i  ][2*j  ]) >= 1.e-15){
                int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL[i];
                int dof_j = 3*numNodesCoarse + 3*numNodesFine + 2*connecL[j];
                MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinA1[2*i  ][2*j  ],ADD_VALUES);
                dof_i++; dof_j++;
                MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinA1[2*i+1][2*j+1],ADD_VALUES);
            };
        };
        //RHS VECTOR
        //COUPLING OPERATOR
        int d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL[i];
        VecSetValues(b,1,&d_i,&Rhs2[2*i  ],ADD_VALUES);

        d_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL[i] + 1;
        VecSetValues(b,1,&d_i,&Rhs2[2*i+1],ADD_VALUES);

        int dof_i = 3*connecC[i];
        VecSetValues(b,1,&dof_i,&rhsLagMult2[2*i  ],ADD_VALUES);

        dof_i = 3*connecC[i] + 1;
        VecSetValues(b,1,&dof_i,&rhsLagMult2[2*i+1],ADD_VALUES);

        //SUPG STABILIZATION
        dof_i = 3*connecC[i];
        VecSetValues(b,1,&dof_i,&localMV_vec[2*i  ],ADD_VALUES);

        dof_i = 3*connecC[i] + 1;
        VecSetValues(b,1,&dof_i,&localMV_vec[2*i+1],ADD_VALUES);

        //PSPG STABILIZATION
        dof_i = 3*connecC[i] + 2;
        VecSetValues(b,1,&dof_i,&localMV_vec[12+i],ADD_VALUES);

        //ARLEQUIN STABILIZATION
        dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL[i];
        VecSetValues(b,1,&dof_i,&RhsArlequin2[2*i  ],ADD_VALUES);
        dof_i++;
        VecSetValues(b,1,&dof_i,&RhsArlequin2[2*i+1],ADD_VALUES);
    };       


    return;
}

//------------------------------------------------------------------------------
//----------------COMPUTE ARLEQUIN COUPLED NAVIER-STOKES PROBLEM----------------
//------------------------------------------------------------------------------
template<>
void Arlequin<2,2>::assembleArlequinSystem(){

    //Coarse mesh
    for (int jel = 0; jel < numElemCoarse; jel++){   
        if (domDecompCoarse.first[jel] == rank) {
            //Compute Element matrix
            double **matrix;
            double rhs[18] = {};
            matrix = new double*[18]();
            for (int i = 0; i < 18; i++) matrix[i] = new double[18]();

            elementsCoarse_[jel] -> getTransientNavierStokes(matrix,rhs);

            setMatVecValuesCoarseModel(matrix,rhs,elementsCoarse_[jel] -> getConnectivity());
            
            for (int i = 0; i < 18; ++i) delete [] matrix[i];
            delete [] matrix;
        };
    };

    //Fine mesh
    for (int jel = 0; jel < numElemFine; jel++){           
        if (domDecompFine.first[jel] == rank) {
            //Compute Element matrix                    
            double **matrix;
            double rhs[18] = {};
            matrix = new double*[18]();
            for (int i = 0; i < 18; i++) matrix[i] = new double[18]();

            elementsFine_[jel] -> getTransientNavierStokes(matrix,rhs);
            
            setMatVecValuesFineModel(matrix,rhs,elementsFine_[jel] -> getConnectivity());

            for (int i = 0; i < 18; ++i) delete [] matrix[i];
            delete [] matrix;
        };                
    };
      

    //Lagrange Multipliers
    for (int l=0; l< numElemGlueZoneFine; l++){
        int jel = elementsGlueZoneFine_[l];
        if (domDecompFine.first[jel] == rank) {
            
            int *connecC;
            int *connec = elementsFine_[jel] -> getConnectivity();
            int *connecL = glueZoneFine_[l] -> getConnectivity();
            //FINE MESH
            //Matrices
            double **Ajac2;
            double **localMV_mat;
            double **ArlequinA1;
            double **ArlequinA2;
            Ajac2 = new double*[18]();
            localMV_mat = new double*[18]();
            ArlequinA1 = new double*[18]();
            ArlequinA2 = new double*[18]();
            for (int i = 0; i < 18; i++){
                Ajac2[i] = new double[18]();
                localMV_mat[i] = new double[18]();
                ArlequinA1[i] = new double[18]();
                ArlequinA2[i] = new double[18]();
            }
            //Vectors
            double rhsLagMult2[18] = {};
            double Rhs2[18] = {};
            double localMV_vec[18] = {};
            double RhsArlequin2[18] = {};

            // FINE MESH
            //Computes element matrix

            elementsFine_[jel] -> getLagrangeMultipliersSameMesh(Ajac2, rhsLagMult2, Rhs2);
            //PSPG and SUPG stabilizations
            elementsFine_[jel] -> getLagrangeMultipliersSUPG_PSPG_SameMesh(localMV_mat,localMV_vec);

            //Arlequin Stabilization
            elementsFine_[jel] -> getLagrangeMultipliersArlequinSameMesh(ArlequinA1, ArlequinA2, RhsArlequin2);
            
            setMatVecValuesLagMultFineFine(Ajac2,localMV_mat,ArlequinA1,ArlequinA2, 
                                           Rhs2,rhsLagMult2,localMV_vec,RhsArlequin2,
                                           elementsFine_[jel] -> getConnectivity(),
                                           glueZoneFine_[l] -> getConnectivity());

            for (int i = 0; i < 18; ++i) {
                delete [] Ajac2[i];
                delete [] localMV_mat[i];
                delete [] ArlequinA1[i];
                delete [] ArlequinA2[i];
            }
            delete [] Ajac2;
            delete [] localMV_mat;
            delete [] ArlequinA1;
            delete [] ArlequinA2;
            
            //COAESE MESH
            //Counts number of coarse mesh intersecting the fine element
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
            //Compute the Lagrange Multiplier element matrix
            for (int ielem = 0; ielem < numElemIntersect; ielem++){
                
                int iElemCoarse = diffElem[ielem];
                double pspg = 0;//elementsCoarse_[iElemCoarse] -> getPSPG();
                double press_[6], velX_[6], velY_[6], velXPrev_[6], velYPrev_[6];

                connecC = elementsCoarse_[iElemCoarse] -> getConnectivity();

                for (int k = 0; k < 6; k++){
                    press_[k] = nodesCoarse_[connecC[k]] -> getPressure();
                    velX_[k] = nodesCoarse_[connecC[k]] -> getVelocity(0);
                    velY_[k] = nodesCoarse_[connecC[k]] -> getVelocity(1);
                    velXPrev_[k] = nodesCoarse_[connecC[k]] -> getPreviousVelocity(0);
                    velYPrev_[k] = nodesCoarse_[connecC[k]] -> getPreviousVelocity(1);
                }
                
                Ajac2 = new double*[18]();
                localMV_mat = new double*[18]();
                ArlequinA1 = new double*[18]();
                ArlequinA2 = new double*[18]();
                for (int i = 0; i < 18; i++){
                    Ajac2[i] = new double[18]();
                    localMV_mat[i] = new double[18]();
                    ArlequinA1[i] = new double[18]();
                    ArlequinA2[i] = new double[18]();
                }
                //Vectors
                double rhsLagMult2[18] = {};
                double Rhs2[18] = {};
                double localMV_vec[18] = {};
                double RhsArlequin2[18] = {};

                elementsFine_[jel] -> getLagrangeMultipliersDifferentMesh(iElemCoarse,pspg,press_,velX_,velY_,velXPrev_,velYPrev_,Ajac2,rhsLagMult2,Rhs2);

                elementsFine_[jel] -> getLagrangeMultipliersSUPG_PSPG_DifferentMesh(iElemCoarse,pspg,press_,velX_,velY_,localMV_mat,localMV_vec);
        
                elementsFine_[jel] -> getLagrangeMultipliersArlequinDifferentMesh(iElemCoarse,pspg,press_,velX_,velY_,ArlequinA1,ArlequinA2,RhsArlequin2);

                setMatVecValuesLagMultFineCoarse(Ajac2, localMV_mat, ArlequinA1, ArlequinA2, 
                                                 Rhs2, rhsLagMult2, localMV_vec, RhsArlequin2,
                                                 elementsCoarse_[iElemCoarse] -> getConnectivity(), 
                                                 glueZoneFine_[l] -> getConnectivity());


                for (int i = 0; i < 18; ++i) {
                    delete [] Ajac2[i];
                    delete [] localMV_mat[i];
                    delete [] ArlequinA1[i];
                    delete [] ArlequinA2[i];
                }
                delete [] Ajac2;
                delete [] localMV_mat;
                delete [] ArlequinA1;
                delete [] ArlequinA2;

            }; //Number of intersections
        }; // if element belongs to the glue zone
    }; // Glue zone

    return;
}


//------------------------------------------------------------------------------
//----------------COMPUTE ARLEQUIN COUPLED NAVIER-STOKES PROBLEM----------------
//------------------------------------------------------------------------------
template<>
void Arlequin<2,2>::initialAcceleration(){

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
            
    std::clock_t t1 = std::clock();
    
    // Preallocates the matrix
    ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                        sysSize, sysSize, 300, NULL, 900, NULL, &A); 
    
    //Create PETSc vectors
    ierr = VecCreate(PETSC_COMM_WORLD, &b); 
    ierr = VecSetSizes(b, PETSC_DECIDE, sysSize); 
    ierr = VecSetFromOptions(b);
    ierr = VecDuplicate(b, &u);
                
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
            double **matrix;
            double rhs[18] = {};
            matrix = new double*[18]();
            for (int i = 0; i < 18; i++) matrix[i] = new double[18]();

            elementsCoarse_[jel] -> getTransientNavierStokesInitial(rhs,matrix);

            setMatVecValuesCoarseModel(matrix,rhs,elementsCoarse_[jel] -> getConnectivity());

            for (int i = 0; i < 18; ++i) delete [] matrix[i];
            delete [] matrix;
        };//Domain decomposition Coarse
    };//Elements Coarse
          
    //Fine mesh
    for (int jel = 0; jel < numElemFine; jel++){   
        
        if (domDecompFine.first[jel] == rank) {
                              
            //Compute Element matrix
            double **matrix;
            double rhs[18] = {};
            matrix = new double*[18]();
            for (int i = 0; i < 18; i++) matrix[i] = new double[18]();

            elementsFine_[jel] -> getTransientNavierStokesInitial(rhs,matrix);
            
            setMatVecValuesFineModel(matrix,rhs,elementsFine_[jel] -> getConnectivity());          
            
            for (int i = 0; i < 18; ++i) delete [] matrix[i];
            delete [] matrix;
        };//Domain decomposition Fine             
    };//Elements Fine
                              
    //Lagrange Multipliers
    for (int l=0; l< numElemGlueZoneFine; l++){
                          
        int jel = elementsGlueZoneFine_[l];
                          
        if (domDecompFine.first[jel] == rank) {              
            int *connecC;
            int *connec = elementsFine_[jel] -> getConnectivity();
            int *connecL = glueZoneFine_[l] -> getConnectivity();

            //FINE MESH
            //Matrices
            double **Ajac2;
            double **localMV_mat;
            double **ArlequinA1;
            double **ArlequinA2;
            Ajac2 = new double*[18]();
            localMV_mat = new double*[18]();
            ArlequinA1 = new double*[18]();
            ArlequinA2 = new double*[18]();
            for (int i = 0; i < 18; i++){
                Ajac2[i] = new double[18]();
                localMV_mat[i] = new double[18]();
                ArlequinA1[i] = new double[18]();
                ArlequinA2[i] = new double[18]();
            }
            //Vectors
            double rhsLagMult2[18] = {};
            double Rhs2[18] = {};
            double localMV_vec[18] = {};
            double RhsArlequin2[18] = {};

            // FINE MESH
            //Computes element matrix

            elementsFine_[jel] -> getLagrangeMultipliersSameMesh(Ajac2, rhsLagMult2, Rhs2);
            //PSPG and SUPG stabilizations
            elementsFine_[jel] -> getLagrangeMultipliersSUPG_PSPG_SameMesh(localMV_mat,localMV_vec);

            //Arlequin Stabilization
            elementsFine_[jel] -> getLagrangeMultipliersArlequinSameMesh(ArlequinA1, ArlequinA2, RhsArlequin2);
            
            setMatVecValuesLagMultFineFine(Ajac2,localMV_mat,ArlequinA1,ArlequinA2, 
                                           Rhs2,rhsLagMult2,localMV_vec,RhsArlequin2,
                                           elementsFine_[jel] -> getConnectivity(),
                                           glueZoneFine_[l] -> getConnectivity());

            for (int i = 0; i < 18; ++i) {
                delete [] Ajac2[i];
                delete [] localMV_mat[i];
                delete [] ArlequinA1[i];
                delete [] ArlequinA2[i];
            }
            delete [] Ajac2;
            delete [] localMV_mat;
            delete [] ArlequinA1;
            delete [] ArlequinA2;

              
            //COAESE MESH
            //Counts number of coarse mesh intersecting the fine element
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
            //Compute the Lagrange Multiplier element matrix
            for (int ielem = 0; ielem < numElemIntersect; ielem++){
                
                int iElemCoarse = diffElem[ielem];
                double pspg = 0;//elementsCoarse_[iElemCoarse] -> getPSPG();
                double press_[6], velX_[6], velY_[6], velXPrev_[6], velYPrev_[6];

                connecC = elementsCoarse_[iElemCoarse] -> getConnectivity();

                for (int k = 0; k < 6; k++){
                    press_[k] = nodesCoarse_[connecC[k]] -> getPressure();
                    velX_[k] = nodesCoarse_[connecC[k]] -> getVelocity(0);
                    velY_[k] = nodesCoarse_[connecC[k]] -> getVelocity(1);
                    velXPrev_[k] = nodesCoarse_[connecC[k]] -> getPreviousVelocity(0);
                    velYPrev_[k] = nodesCoarse_[connecC[k]] -> getPreviousVelocity(1);
                }
                
                Ajac2 = new double*[18]();
                localMV_mat = new double*[18]();
                ArlequinA1 = new double*[18]();
                ArlequinA2 = new double*[18]();
                for (int i = 0; i < 18; i++){
                    Ajac2[i] = new double[18]();
                    localMV_mat[i] = new double[18]();
                    ArlequinA1[i] = new double[18]();
                    ArlequinA2[i] = new double[18]();
                }
                //Vectors
                double rhsLagMult2[18] = {};
                double Rhs2[18] = {};
                double localMV_vec[18] = {};
                double RhsArlequin2[18] = {};

                elementsFine_[jel] -> getLagrangeMultipliersDifferentMesh(iElemCoarse,pspg,press_,velX_,velY_,velXPrev_,velYPrev_,Ajac2,rhsLagMult2,Rhs2);

                elementsFine_[jel] -> getLagrangeMultipliersSUPG_PSPG_DifferentMesh(iElemCoarse,pspg,press_,velX_,velY_,localMV_mat,localMV_vec);
        
                elementsFine_[jel] -> getLagrangeMultipliersArlequinDifferentMesh(iElemCoarse,pspg,press_,velX_,velY_,ArlequinA1,ArlequinA2,RhsArlequin2);

                setMatVecValuesLagMultFineCoarse(Ajac2, localMV_mat, ArlequinA1, ArlequinA2, 
                                                 Rhs2, rhsLagMult2, localMV_vec, RhsArlequin2,
                                                 elementsCoarse_[iElemCoarse] -> getConnectivity(), 
                                                 glueZoneFine_[l] -> getConnectivity());


                for (int i = 0; i < 18; ++i) {
                    delete [] Ajac2[i];
                    delete [] localMV_mat[i];
                    delete [] ArlequinA1[i];
                    delete [] ArlequinA2[i];
                }
                delete [] Ajac2;
                delete [] localMV_mat;
                delete [] ArlequinA1;
                delete [] ArlequinA2;

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
        // if (nodesCoarse_[i] -> getDistFunction() > -1.2) val *= 1000.e0;
        u_[0] = nodesCoarse_[i] -> getVelocity(0);
        normU += val*val;
        nodesCoarse_[i] -> setAccelerationComponent(0,val);
        nodesCoarse_[i] -> setPreviousVelocityComponent(0,u_[0]-val*gamma*dTime);
          
        Ii = 2 * i + 1;
        ierr = VecGetValues(All, Ione, &Ii, &val);
        // if (nodesCoarse_[i] -> getDistFunction() > -1.2) val *= 1000.e0;
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
                                        
    std::clock_t t2 = std::clock();
      
    ierr = KSPDestroy(&ksp); 
    ierr = VecDestroy(&b); 
    ierr = VecDestroy(&u); 
    ierr = VecDestroy(&All); 
    ierr = MatDestroy(&A); 
    // ierr = MatDestroy(&F); CHKERRQ(ierr);
    //ierr = MatDestroy(&C); CHKERRQ(ierr);
    
    //Printing results
    if (rank == 0) {
        printResultsCoarse(100);
        printResultsFine(100);
    }
        
    return;

};






//------------------------------------------------------------------------------
//----------------COMPUTE ARLEQUIN COUPLED NAVIER-STOKES PROBLEM----------------
//------------------------------------------------------------------------------
template<>
int Arlequin<2,2>::solveArlequinProblem(int iterNumber, double tolerance,
                                      int problem_type, int time_dependency){


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
    if (rank == 0) {
        printResultsCoarse(100);
        printResultsFine(100);
    }
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
        // std::cout << " asd " << iTimeStep << " " << rank << std::endl;

        if(iTimeStep == 10){
            double dd = 0.75;
            parametersFine -> setSpectralRadius(dd);
            parametersCoarse -> setSpectralRadius(dd);
            std::cout << "AQUI " << rank << std::endl;
        }
           
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
            
            std::clock_t t1 = std::clock();
            
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
                        
            // for (int i=0; i<sysSize; i++){
            //     double val = 1.e-20;
            //     ierr = MatSetValues(A,1,&i,1,&i,&val,ADD_VALUES);
                
            // }
            

            std::clock_t t3 = std::clock();
            assembleArlequinSystem();
            std::clock_t t4 = std::clock();

            //std::cout << "Enter PETSc " << rank << std::endl;
            
            //Assemble matrices and vectors
            ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
            ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
            
            ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
            ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
            std::clock_t t5 = std::clock();

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
            // ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);



            
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
            // std::cout << "AQQQEQE1 " << rank << std::endl;
            // ierr = PCSetFromOptions(pc);CHKERRQ(ierr);
            // std::cout << "AQQQEQE2 " << rank << std::endl;
            ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
            // std::cout << "AQQQEQE3 " << rank << std::endl;
            ierr = KSPSetUp(ksp);
            // std::cout << "AQQQEQE4 " << rank << std::endl;
            
// #if defined(PETSC_HAVE_MUMPS)
//             PetscInt  info1,info2,icntl14;
            
//             MatMumpsGetInfo(F,1,&info1);
//             MatMumpsGetInfo(F,2,&info2);

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
//             MatMumpsGetIcntl(F,14,&icntl14);    
//             if((rank == 0) && (info1 != 0)) std::cout << " INFO(1) = " << info1
//                                                       << " " << info2 << " " 
//                                                       << icntl14 << std::endl;
// #endif
            
            // ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
            
            ierr = KSPSolve(ksp,b,u);CHKERRQ(ierr);
            
            ierr = KSPGetTotalIterations(ksp, &iterations);
            
            std::clock_t t6 = std::clock();

            if (rank == 0) std::cout << "TIME " << 1000.*(t4-t3)/CLOCKS_PER_SEC/1000. << " " 
                                                << 1000.*(t5-t4)/CLOCKS_PER_SEC/1000. << " " 
                                                << 1000.*(t6-t5)/CLOCKS_PER_SEC/1000. << std::endl;
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
                Ii = 3 * i;
                ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                // if (nodesCoarse_[i] -> getDistFunction() > -1.2) val *= 1000.e0;
                u_[0] = val;
                normU += val*val;
                nodesCoarse_[i] -> incrementAcceleration(0,u_[0]);
                nodesCoarse_[i] -> incrementVelocity(0,u_[0]*gamma*dTime);
                
                Ii = 3 * i + 1;
                ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                // if (nodesCoarse_[i] -> getDistFunction() > -1.2) val *= 1000.e0;
                u_[1] = val;
                normU += val*val;
                nodesCoarse_[i] -> incrementAcceleration(1,u_[1]);
                nodesCoarse_[i] -> incrementVelocity(1,u_[1]*gamma*dTime);

                Ii = 3 * i + 2;
                ierr = VecGetValues(All,Ione,&Ii,&val);CHKERRQ(ierr);
                p_ = val;
                normP += val*val;
                nodesCoarse_[i] -> incrementPressure(p_);
            };
            
            for (int i = 0; i < numNodesFine; ++i){
                Ii = 3 * numNodesCoarse + 3 * i;
                ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                u_[0] = val;
                normU += val*val;
                nodesFine_[i] -> incrementAcceleration(0,u_[0]);
                nodesFine_[i] -> incrementVelocity(0,u_[0]*gamma*dTime);
                
                Ii = 3 * numNodesCoarse + 3 * i + 1;
                ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                u_[1] = val;
                normU += val*val;
                nodesFine_[i] -> incrementAcceleration(1,u_[1]);
                nodesFine_[i] -> incrementVelocity(1,u_[1]*gamma*dTime);
        
                Ii = 3 * numNodesCoarse + 3 * i + 2;
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
            
            std::clock_t t2 = std::clock();
            
            if(rank == 0){
                std::cout<<"Iteration = " << inewton << " (" << iterations <<  
                    ")  Du Norm = " << std::scientific << sqrt(normU) 
                         << " " << sqrt(normP) 
                         << " " << sqrt(normL) 
                         << " " << val << 
                    "  Time (s) = " << std::fixed << 
                    1000.*(t2-t1)/CLOCKS_PER_SEC/1000. << std::endl;
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
        QuadShapeFunction<2,2>                       shapeQuad;
        double phi_[6] = {};
        
        for (int i = 0; i<numNodesFine; i++){
            nodesFine_[i] -> setVelocityArlequin(0,nodesFine_[i] -> 
                                                 getVelocity(0));
            nodesFine_[i] -> setVelocityArlequin(1,nodesFine_[i] -> 
                                                 getVelocity(1));
            nodesFine_[i] -> setPressureArlequin(nodesFine_[i] ->getPressure());
        };
        
        for (int i = 0; i<numNodesGlueZoneFine; i++){
            double u_coarse[6], v_coarse[6], p_coarse[6];
            int *connecCoarse;
            
            double u = 0.;
            double v = 0.;
            double p = 0.;
            
            int elCoarse = nodesFine_[nodesGlueZoneFine_[i]] -> getNodalElemCorrespondence();
            double* xsi = nodesFine_[nodesGlueZoneFine_[i]] -> getNodalXsiCorrespondence();
                        
            connecCoarse = elementsCoarse_[elCoarse] -> getConnectivity();
            
            for (int j=0; j<6; j++){
                u_coarse[j] = nodesCoarse_[connecCoarse[j]] -> getVelocity(0);
                v_coarse[j] = nodesCoarse_[connecCoarse[j]] -> getVelocity(1);
                p_coarse[j] = nodesCoarse_[connecCoarse[j]] -> getPressure();
            };
            
            shapeQuad.evaluate(xsi,phi_);
            
            for (int j=0; j<6; j++){
                u += u_coarse[j] * phi_[j];
                v += v_coarse[j] * phi_[j];
                p += p_coarse[j] * phi_[j];
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
            //Printing results
            printResultsCoarse(iTimeStep);
            printResultsFine(iTimeStep);
        };        
 

    };
        
    return 0;

};



//------------------------------------------------------------------------------
//----------------COMPUTE ARLEQUIN COUPLED NAVIER-STOKES PROBLEM----------------
//------------------------------------------------------------------------------
template<>
int Arlequin<2,2>::solveArlequinProblemMoving(int iterNumber, double tolerance,
                                            int problem_type, 
                                            int time_dependency){


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
    
    numTimeSteps = 2000;

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



            double xn[2];
            double *xi = nodesFine_[i] -> getInitialCoordinates();       
            double *x = nodesFine_[i] -> getCoordinates();       
    
            double a = -20 * pi / 180 + 10 * pi / 180 * std::cos(2.*pi*iTimeStep*dTime);// + 10 * pi / 180;

            // std::cout << " AAA " << a << std::endl;

            xn[0] = 0.5 + (xi[0]-0.5) * std::cos(a) - (xi[1]-0.0) * std::sin(a);
            xn[1] = 0.0 + (xi[0]-0.5) * std::sin(a) + (xi[1]-0.0) * std::cos(a);

            u[0] = (xn[0] - x[0]) / dTime;
            u[1] = (xn[1] - x[1]) / dTime;


            nodesFine_[i] -> setMeshVelocity(u);
      
            nodesFine_[i] -> setPreviousCoordinates(0,x[0]);
            nodesFine_[i] -> setPreviousCoordinates(1,x[1]);

            nodesFine_[i] -> setCoordinates(xn);
        };
        
        setNodalCorrespondenceFine();
        setSignaledDistance();
        setWeightFunction(1.);
        


        //STARTS NEWTON-RAPHSON
        for (int inewton = 0; inewton < iterNumber; inewton++){
            
            std::clock_t t1 = std::clock();
            
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
                        
            for (int i=0; i<sysSize; i++){
                double val = 1.e-15;
                ierr = MatSetValues(A,1,&i,1,&i,&val,ADD_VALUES);
                
            }

            assembleArlequinSystem();
            
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
                Ii = 3 * i;
                ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                // if (nodesCoarse_[i] -> getDistFunction() > -1.2) val *= 1000.e0;
                u_[0] = val;
                normU += val*val;
                nodesCoarse_[i] -> incrementAcceleration(0,u_[0]);
                nodesCoarse_[i] -> incrementVelocity(0,u_[0]*gamma*dTime);
                
                Ii = 3 * i + 1;
                ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                // if (nodesCoarse_[i] -> getDistFunction() > -1.2) val *= 1000.e0;
                u_[1] = val;
                normU += val*val;
                nodesCoarse_[i] -> incrementAcceleration(1,u_[1]);
                nodesCoarse_[i] -> incrementVelocity(1,u_[1]*gamma*dTime);

                Ii = 3 * i + 2;
                ierr = VecGetValues(All,Ione,&Ii,&val);CHKERRQ(ierr);
                p_ = val;
                normP += val*val;
                nodesCoarse_[i] -> incrementPressure(p_);
            };
            
            for (int i = 0; i < numNodesFine; ++i){
                Ii = 3 * numNodesCoarse + 3 * i;
                ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                u_[0] = val;
                normU += val*val;
                nodesFine_[i] -> incrementAcceleration(0,u_[0]);
                nodesFine_[i] -> incrementVelocity(0,u_[0]*gamma*dTime);
                
                Ii = 3 * numNodesCoarse + 3 * i + 1;
                ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                u_[1] = val;
                normU += val*val;
                nodesFine_[i] -> incrementAcceleration(1,u_[1]);
                nodesFine_[i] -> incrementVelocity(1,u_[1]*gamma*dTime);
        
                Ii = 3 * numNodesCoarse + 3 * i + 2;
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
            
            std::clock_t t2 = std::clock();
            
            if(rank == 0){                
                std::cout<<"Iteration = " << inewton << " (" << iterations <<  
                    ")  Du Norm = " << std::scientific << sqrt(normU) 
                         << " " << sqrt(normP) 
                         << " " << sqrt(normL) 
                         << " " << val << 
                    "  Time (s) = " << std::fixed << 
                    1000.*(t2-t1)/CLOCKS_PER_SEC/1000. << std::endl;
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
        
        QuadShapeFunction<2,2>                       shapeQuad;
        double phi_[6] = {};
        
        for (int i = 0; i<numNodesFine; i++){
            nodesFine_[i] -> setVelocityArlequin(0,nodesFine_[i] -> 
                                                 getVelocity(0));
            nodesFine_[i] -> setVelocityArlequin(1,nodesFine_[i] -> 
                                                 getVelocity(1));
            nodesFine_[i] -> setPressureArlequin(nodesFine_[i] ->getPressure());
        };
        
        for (int i = 0; i<numNodesGlueZoneFine; i++){
            double u_coarse[6], v_coarse[6], p_coarse[6];
            
            double u = 0.;
            double v = 0.;
            double p = 0.;
            
            int elCoarse = nodesFine_[nodesGlueZoneFine_[i]] -> getNodalElemCorrespondence();
            double* xsi = nodesFine_[nodesGlueZoneFine_[i]] -> getNodalXsiCorrespondence();
            
            int *connecCoarse = elementsCoarse_[elCoarse] -> getConnectivity();
            
            for (int j=0; j<6; j++){
                u_coarse[j] = nodesCoarse_[connecCoarse[j]] -> getVelocity(0);
                v_coarse[j] = nodesCoarse_[connecCoarse[j]] -> getVelocity(1);
                p_coarse[j] = nodesCoarse_[connecCoarse[j]] -> getPressure();
            };
            
            shapeQuad.evaluate(xsi,phi_);
            
            for (int j=0; j<6; j++){
                u += u_coarse[j] * phi_[j];
                v += v_coarse[j] * phi_[j];
                p += p_coarse[j] * phi_[j];
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
                       
            //Printing results
            printResultsCoarse(iTimeStep);
            printResultsFine(iTimeStep);
        };
    };
        
    return 0;

};



//------------------------------------------------------------------------------
//----------------COMPUTE ARLEQUIN COUPLED NAVIER-STOKES PROBLEM----------------
//------------------------------------------------------------------------------
template<>
int Arlequin<2,2>::solveFSIArlequin(int iterNumber, double tolerance,
                                  int problem_type, int iTimeStep){


    std::ofstream dragLift;
    dragLift.open("dragLift.dat", std::ofstream::out | std::ofstream::app);

    if ((problem_type < 1) || (problem_type > 2)){
        std::cout << "WRONG PROBLEM TYPE." << std::endl;
        return 0;
    };

    // Computes the system size
    int sysSize = 3 * numNodesCoarse +  
        + 3 * numNodesFine + 2 * numNodesGlueZoneFine;
    
    setNodalCorrespondenceFine();
    setSignaledDistance();
    setWeightFunction(1.);

    // initialAcceleration();

    double &alpha_f = parametersFine -> getAlphaF();
    double &alpha_m = parametersFine -> getAlphaM();
    double &gamma = parametersFine -> getGamma();
    
    //STARTS NEWTON-RAPHSON
    for (int inewton = 0; inewton < iterNumber; inewton++){
        
        std::clock_t t1 = std::clock();
        
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
                    
        for (int i=0; i<sysSize; i++){
            double val = 1.e-20;
            ierr = MatSetValues(A,1,&i,1,&i,&val,ADD_VALUES);
            
        }
                            
        assembleArlequinSystem();
        
        //Assemble matrices and vectors
        ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        
        ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
        ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
        
        //Create KSP context to solve the linear system
        ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
        
        ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
        



        // // ierr = MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE,0, NULL, &nullsp);
        // // ierr = MatSetNullSpace(A, nullsp);
        // // ierr = MatNullSpaceDestroy(&nullsp);


        // // if (iTimeStep > 5){

        //     ierr = KSPSetTolerances(ksp,1.e-7,1.e-10,PETSC_DEFAULT,
        //                             15);CHKERRQ(ierr);
            
        //     //ierr = KSPGMRESSetRestart(ksp, 10); CHKERRQ(ierr);
            
        //     ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
            
        //     ierr = PCSetType(pc,PCJACOBI);CHKERRQ(ierr);
            
        //     //ierr = KSPSetPCSide(ksp, PC_RIGHT);
        //     ierr = KSPSetType(ksp,KSPGMRES); CHKERRQ(ierr);
            
        //     ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
        //     // ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
        // // }else{


        
#if defined(PETSC_HAVE_MUMPS)
        ierr = KSPSetType(ksp,KSPPREONLY);
        ierr = KSPGetPC(ksp,&pc);
        ierr = PCSetType(pc, PCLU);
        
        // ierr = PCFactorSetMatSolverType(pc,MATSOLVERMUMPS);
        // PCFactorSetUpMatSolverType(pc);
        // PCFactorGetMatrix(pc,&F);
        
        // PetscInt ival,icntl;
        // icntl = 14; ival = 80;
        // MatMumpsSetIcntl(F,icntl,ival);
        // icntl = 28; ival = 2;
        // MatMumpsSetIcntl(F,icntl,ival);
        // icntl = 29; ival = 2;
        // MatMumpsSetIcntl(F,icntl,ival);
        // icntl = 16; ival = 0;
        // MatMumpsSetIcntl(F,icntl,ival);
        // icntl = 4; ival = 3;
        // MatMumpsSetIcntl(F,icntl,ival);
        // icntl = 11; ival = 1;
        // MatMumpsSetIcntl(F,11,1);

        //MatMumpsSetIcntl(F,21,0);

        
#endif
        ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
        ierr = KSPSetUp(ksp);
        
// #if defined(PETSC_HAVE_MUMPS)
//         PetscInt  info1,info2,icntl14;
        
//         MatMumpsGetInfo(F,1,&info1);
//         MatMumpsGetInfo(F,2,&info2);

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
      
        // // if(rank==0) std::cout << "INFOG = " << info5 << " " << info6 << " " << info7 << " " << info8 << " " << info9 << " " << info10 << " " << info11 << std::endl;
        // MatMumpsGetIcntl(F,14,&icntl14);    
        // if((rank == 0) && (info1 != 0)) std::cout << " INFO(1) = " << info1
        //                                           << " " << info2 << " " 
        //                                           << icntl14 << std::endl;
// #endif
        // }
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
            Ii = 3 * i;
            ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
            // if (nodesCoarse_[i] -> getDistFunction() > -1.2) val *= 1000.e0;
            u_[0] = val;
            normU += val*val;
            nodesCoarse_[i] -> incrementAcceleration(0,u_[0]);
            nodesCoarse_[i] -> incrementVelocity(0,u_[0]*gamma*dTime);
            
            Ii = 3 * i + 1;
            ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
            // if (nodesCoarse_[i] -> getDistFunction() > -1.2) val *= 1000.e0;
            u_[1] = val;
            normU += val*val;
            nodesCoarse_[i] -> incrementAcceleration(1,u_[1]);
            nodesCoarse_[i] -> incrementVelocity(1,u_[1]*gamma*dTime);

            Ii = 3 * i + 2;
            ierr = VecGetValues(All,Ione,&Ii,&val);CHKERRQ(ierr);
            p_ = val;
            normP += val*val;
            nodesCoarse_[i] -> incrementPressure(p_);
        };
        
        for (int i = 0; i < numNodesFine; ++i){
            Ii = 3 * numNodesCoarse + 3 * i;
            ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
            u_[0] = val;
            normU += val*val;
            nodesFine_[i] -> incrementAcceleration(0,u_[0]);
            nodesFine_[i] -> incrementVelocity(0,u_[0]*gamma*dTime);
            
            Ii = 3 * numNodesCoarse + 3 * i + 1;
            ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
            u_[1] = val;
            normU += val*val;
            nodesFine_[i] -> incrementAcceleration(1,u_[1]);
            nodesFine_[i] -> incrementVelocity(1,u_[1]*gamma*dTime);
    
            Ii = 3 * numNodesCoarse + 3 * i + 2;
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
        
        std::clock_t t2 = std::clock();
        
        if(rank == 0){
            std::cout<<"Iteration = " << inewton << " (" << iterations <<  
                ")  Du Norm = " << std::scientific << sqrt(normU) 
                     << " " << sqrt(normP) 
                     << " " << sqrt(normL) 
                     << " " << val << 
                "  Time (s) = " << std::fixed << 
                1000.*(t2-t1)/CLOCKS_PER_SEC/1000. << std::endl;
        };
        
        ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
        ierr = VecDestroy(&b); CHKERRQ(ierr);
        ierr = VecDestroy(&u); CHKERRQ(ierr);
        ierr = VecDestroy(&All); CHKERRQ(ierr);
        ierr = MatDestroy(&A); CHKERRQ(ierr);
        // ierr = MatDestroy(&F); CHKERRQ(ierr);
        
        
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


    // std::cout << "AQUI1 " << rank << std::endl;


    //Compute real velocity
    QuadShapeFunction<2,2>                       shapeQuad;
    double phi_[6] = {};
    
    for (int i = 0; i<numNodesFine; i++){
        nodesFine_[i] -> setVelocityArlequin(0,nodesFine_[i] -> 
                                             getVelocity(0));
        nodesFine_[i] -> setVelocityArlequin(1,nodesFine_[i] -> 
                                             getVelocity(1));
        nodesFine_[i] -> setPressureArlequin(nodesFine_[i] ->getPressure());
    };
    
    for (int i = 0; i<numNodesGlueZoneFine; i++){
        double u_coarse[6], v_coarse[6], p_coarse[6];
        int *connecCoarse;
        
        double u = 0.;
        double v = 0.;
        double p = 0.;
        
        int elCoarse = nodesFine_[nodesGlueZoneFine_[i]] -> getNodalElemCorrespondence();
        double* xsi = nodesFine_[nodesGlueZoneFine_[i]] -> getNodalXsiCorrespondence();
        
        connecCoarse = elementsCoarse_[elCoarse] -> getConnectivity();
        
        for (int j=0; j<6; j++){
            u_coarse[j] = nodesCoarse_[connecCoarse[j]] -> getVelocity(0);
            v_coarse[j] = nodesCoarse_[connecCoarse[j]] -> getVelocity(1);
            p_coarse[j] = nodesCoarse_[connecCoarse[j]] -> getPressure();
        };
        
        shapeQuad.evaluate(xsi,phi_);
        
        for (int j=0; j<6; j++){
            u += u_coarse[j] * phi_[j];
            v += v_coarse[j] * phi_[j];
            p += p_coarse[j] * phi_[j];
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
    // std::cout << "AQUI2 " << rank << std::endl;

    for (int i=0; i<numNodesCoarse; i++){
        nodesCoarse_[i] -> 
            setVelocityArlequin(0,nodesCoarse_[i] -> getVelocity(0));
        nodesCoarse_[i] -> 
            setVelocityArlequin(1,nodesCoarse_[i] -> getVelocity(1));
        nodesCoarse_[i] -> 
            setPressureArlequin(nodesCoarse_[i] -> getPressure());
    };
    // std::cout << "AQUI3 " << rank << std::endl;
    
    return 0;

};


#endif
