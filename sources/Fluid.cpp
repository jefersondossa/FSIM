#include "Fluid.h"
#include "petscpartitioner.h"

//------------------------------------------------------------------------------
//----------------------COMPUTES DRAG AND LIFT COEFFICIENTS---------------------
//------------------------------------------------------------------------------
void Fluid::dragAndLiftCoefficients(std::ofstream& dragLift){

    double dragCoefficient = 0.;
    double liftCoefficient = 0.;
    double pressureDragCoefficient = 0.;
    double pressureLiftCoefficient = 0.;
    double frictionDragCoefficient = 0.;
    double frictionLiftCoefficient = 0.;
    
    for (int jel = 0; jel < NBoundElements(); jel++){   
        
        double rhoInf = 1.0;
        double velocityInf[2];
        velocityInf[0] = -1.;
        velocityInf[1] = 0.;
        
        double dForce = 0.;
        double lForce = 0.;
        double pDForce = 0.;
        double pLForce = 0.;
        double fDForce = 0.;
        double fLForce = 0.;
        double aux_Mom = 0.;
        double aux_Per = 0.;
        
       for (int i=0; i<numberOfLines; i++){
            //std::cout << "Bound group " << BoundaryVec()[jel] -> getBoundaryGroup() << std::endl;
            if (BoundaryVec()[jel] -> getBoundaryGroup() == dragAndLiftBoundary[i]){
                //std::cout << "AQUI " << numberOfLines<< " " << i << " " << dragAndLiftBoundary[i] << std::endl;
                int iel = BoundaryVec()[jel] -> getElement();
                // ElementVec()[iel] -> computeDragAndLiftForces(pDForce, pLForce, fDForce, fLForce, dForce, lForce, aux_Mom, aux_Per);
                // ElementVec()[iel] -> computeSeparationAngle();
            };
        };
        
        pressureDragCoefficient += pDForce ;/// 
            //(0.5 * rhoInf * velocityInf[0] * velocityInf[0]);
        pressureLiftCoefficient += pLForce;// / 
            //(0.5 * rhoInf * velocityInf[0] * velocityInf[0]);
        
        frictionDragCoefficient += fDForce / 
            (0.5 * rhoInf * velocityInf[0] * velocityInf[0]);
        frictionLiftCoefficient += fLForce / 
            (0.5 * rhoInf * velocityInf[0] * velocityInf[0]);
        
        dragCoefficient += dForce / 
            (0.5 * rhoInf * velocityInf[0] * velocityInf[0]);
        liftCoefficient += lForce / 
            (0.5 * rhoInf * velocityInf[0] * velocityInf[0]);
        
    };
    // std::cout << "vazao " << pressureDragCoefficient << " " << pressureLiftCoefficient << std::endl;
    if (rank == 0) {
        const int timeWidth = 11;
        const int numWidth = 11;
        dragLift << std::setprecision(3) << std::scientific;
        dragLift << std::left << std::setw(timeWidth) ;
        dragLift << std::setw(numWidth) << pressureDragCoefficient;
        dragLift << std::setw(numWidth) << pressureLiftCoefficient;
        dragLift << std::setw(numWidth) << frictionDragCoefficient;
        dragLift << std::setw(numWidth) << frictionLiftCoefficient;
        dragLift << std::setw(numWidth) << dragCoefficient;
        dragLift << std::setw(numWidth) << liftCoefficient;
        dragLift << std::endl;
    }
}

//------------------------------------------------------------------------------
//-------------------------SOLVE STEADY LAPLACE PROBLEM-------------------------
//------------------------------------------------------------------------------

int Fluid::solveSteadyLaplaceProblem(int iterNumber, double tolerance) {

    Vec               b, u, All;
    PetscErrorCode    ierr;
    PetscInt          Istart, Iend, Ii, Ione, iterations;
    KSP               ksp;
    PC                pc;
    VecScatter        ctx;
    PetscScalar       val;
   
    int rank;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);        

    for (int inewton = 0; inewton < iterNumber; inewton++){

        ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                            2*NNodes(), 2*NNodes(),
                            100,NULL,100,NULL,&A); 
        
        ierr = MatGetOwnershipRange(A, &Istart, &Iend);
        
        //Create PETSc vectors
        ierr = VecCreate(PETSC_COMM_WORLD,&b);
        ierr = VecSetSizes(b,PETSC_DECIDE,2*NNodes());
        ierr = VecSetFromOptions(b);
        ierr = VecDuplicate(b,&u);
        ierr = VecDuplicate(b,&All);
        
        //std::cout << "Istart = " << Istart << " Iend = " << Iend << std::endl;
        
        for (int jel = 0; jel < NElements(); jel++){               
            //Compute Element matrix
            VecInt connec = ElementVec()[jel] -> getConnectivity();

            MatrixDouble matrix(NLocDOF(),NLocDOF());
            VecDouble rhs(NLocDOF());
            rhs.setZero();
            matrix.setZero();

            ElementVec()[jel] -> getSolidProblem(matrix,rhs);
            
            //Disperse local contributions into the global matrix
            //Matrix K and C
            for (int i=0; i<NElNodes(); i++){
                for (int j=0; j<NElNodes(); j++){
                    if (fabs(matrix(2*i  ,2*j  )) >= 1.e-8){
                        int dof_i = 2*connec[i];
                        int dof_j = 2*connec[j];
                        ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&matrix(2*i  ,2*j  ),ADD_VALUES);
                    };
                    if (fabs(matrix(2*i+1,2*j  )) >= 1.e-8){
                        int dof_i = 2*connec[i]+1;
                        int dof_j = 2*connec[j];
                        ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&matrix(2*i+1,2*j  ),ADD_VALUES);
                    };
                    if (fabs(matrix(2*i  ,2*j+1)) >= 1.e-8){
                        int dof_i = 2*connec[i];
                        int dof_j = 2*connec[j]+1;
                        ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&matrix(2*i  ,2*j+1),ADD_VALUES);
                    };
                    if (fabs(matrix(2*i+1,2*j+1)) >= 1.e-8){
                        int dof_i = 2*connec[i]+1;
                        int dof_j = 2*connec[j]+1;
                        ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&matrix(2*i+1,2*j+1),ADD_VALUES);
                    };
                };
                                    
                //Rhs vector
                if (fabs(rhs[2*i  ]) >= 1.e-8){
                    int dof_i = 2*connec[i];
                    ierr = VecSetValues(b,1,&dof_i,&rhs[2*i  ],ADD_VALUES);
                };
                
                if (fabs(rhs[2*i+1]) >= 1.e-8){
                    int dof_i = 2*connec[i]+1;
                    ierr = VecSetValues(b,1,&dof_i,&rhs[2*i+1],ADD_VALUES);
                };
            };
        };
        
        //Assemble matrices and vectors
        ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
        ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
        
        ierr = VecAssemblyBegin(b);
        ierr = VecAssemblyEnd(b);
        
        // MatView(A,PETSC_VIEWER_STDOUT_WORLD);
        // VecView(b,PETSC_VIEWER_STDOUT_WORLD);
        
        //Create KSP context to solve the linear system
        ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);
        
        ierr = KSPSetOperators(ksp,A,A);
        
#if defined(PETSC_HAVE_MUMPS)
        ierr = KSPSetType(ksp,KSPPREONLY);
        ierr = KSPGetPC(ksp,&pc);
        ierr = PCSetType(pc, PCLU);
#endif
        
        ierr = KSPSetFromOptions(ksp);
        ierr = KSPSetUp(ksp);
        
        
        
        ierr = KSPSolve(ksp,b,u);
        
        ierr = KSPGetTotalIterations(ksp, &iterations);

        //std::cout << "GMRES Iterations = " << iterations << std::endl;
        
        //ierr = VecView(u,PETSC_VIEWER_STDOUT_WORLD);
        
        //Gathers the solution vector to the master process
        ierr = VecScatterCreateToAll(u, &ctx, &All);
        
        ierr = VecScatterBegin(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);
        
        ierr = VecScatterEnd(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);
        
        ierr = VecScatterDestroy(&ctx);
                
        //Updates nodal values
        double u_ [2];
        Ione = 1;

        for (int i = 0; i < NNodes(); ++i){
            Ii = 2*i;
            ierr = VecGetValues(All, Ione, &Ii, &val);
            u_[0] = val;
            Ii = 2*i+1;
            ierr = VecGetValues(All, Ione, &Ii, &val);
            u_[1] = val;
            if (NodeVec()[i] -> getConstrainsLaplace(0) != 1) NodeVec()[i] -> incrementCoordinate(0,u_[0]);
            if (NodeVec()[i] -> getConstrainsLaplace(1) != 1) NodeVec()[i] -> incrementCoordinate(1,u_[1]);




            // double x = NodeVec()[i] -> getCoordinateValue(0);
            // double y = NodeVec()[i] -> getCoordinateValue(1);
            // double xp = NodeVec()[i] -> getPreviousCoordinateValue(0);
            // double yp = NodeVec()[i] -> getPreviousCoordinateValue(1);
            // double vx = NodeVec()[i] -> getMeshVelocity(0);
            // double vy = NodeVec()[i] -> getMeshVelocity(1);
            // double ax = NodeVec()[i] -> getMeshAcceleration(0);
            // double ay = NodeVec()[i] -> getMeshAcceleration(1);

            // double accelx = (x - xp) / (0.25 * dTime * dTime) - vx / (0.25 * dTime) - ax * (0.5/0.25 - 1.0);
            // double accely = (y - yp) / (0.25 * dTime * dTime) - vy / (0.25 * dTime) - ay * (0.5/0.25 - 1.0);

            // NodeVec()[i] -> setMeshAccelerationComponent(0,accelx);
            // NodeVec()[i] -> setMeshAccelerationComponent(1,accely);

            // double velx = 0.5 * dTime * accelx + vx + dTime * (1.0 - 0.5) * ax;
            // double vely = 0.5 * dTime * accely + vy + dTime * (1.0 - 0.5) * ay;
            // // std::cout << "asd as " << vx << " " << velx << " " << accelx << " " << ax << std::endl;
            // NodeVec()[i] -> setMeshVelocityComponent(0,velx);
            // NodeVec()[i] -> setMeshVelocityComponent(1,vely);

        };
        
        //Computes the solution vector norm
        ierr = VecNorm(u,NORM_2,&val);
        
        if(rank == 0){
            std::cout << "MESH MOVING - ERROR = " << val 
                      << std::scientific <<  std::endl;
        };

        ierr = KSPDestroy(&ksp); 
        ierr = VecDestroy(&b); 
        ierr = VecDestroy(&u); 
        ierr = VecDestroy(&All); 
        ierr = MatDestroy(&A); 

        if(val <= tolerance){
            break;
        };
    };
    
    // for (int i=0; i<NElements(); i++){
    //     ElementVec()[i] -> computeNodalGradient();            
    // };

    if (rank == 0) {
        //Computing velocity divergent
        //      printResults(1);
    };

    return 0;
};


//------------------------------------------------------------------------------
//-------------------------SOLVE TRANSIENT FLUID PROBLEM------------------------
//------------------------------------------------------------------------------

int Fluid::solveFSIFluid(int iterNumber, double tolerance, int problem_type){

    Vec               b, u, All, Allu;
    PetscErrorCode    ierr;
    PetscInt          Ii, Ione, iterations;
    KSP               ksp;
    PC                pc;
    VecScatter        ctx;
    PetscScalar       val;
    //    MatNullSpace      nullsp;

    int rank;

    iAux++;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    //Check if the problem type can be computed
    if ((problem_type > 2) || (problem_type < 1)){
        std::cout << "WRONG PROBLEM TYPE." << std::endl;
        return 0;
    };
            
    double duNorm=100.;

    double &alpha_f = fProbParameters.getAlphaF();
    double &alpha_m = fProbParameters.getAlphaM();
    double &gamma = fProbParameters.getGamma();
         
    for (int inewton = 0; inewton < iterNumber; inewton++){
        
        std::clock_t t1 = std::clock();
        
        ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                            3*NNodes(), 3*NNodes(),
                            100,NULL,300,NULL,&A); 
        
                
        //Create PETSc vectors
        ierr = VecCreate(PETSC_COMM_WORLD,&b);
        ierr = VecSetSizes(b,PETSC_DECIDE,3*NNodes());
        
        ierr = VecSetFromOptions(b);
        ierr = VecDuplicate(b,&u);
        ierr = VecDuplicate(b,&All);
        
        //std::cout << "Istart = " << Istart << " Iend = " << Iend << std::endl;

        for (int jel = 0; jel < NElements(); jel++){   
            
            if (part_elem[jel] == rank) {
                //Compute Element matrix
                VecInt connec = ElementVec()[jel] -> getConnectivity();

                MatrixDouble matrix(18,18);
                matrix.setZero();
                VecDouble rhs(18);
                rhs.setZero();

                ElementVec()[jel] -> getTransientNavierStokes(matrix,rhs);
                
                //Disperse local contributions into the global matrix
                //Matrix K and C
                for (int i=0; i<6; i++){
                    for (int j=0; j<6; j++){
                        int dof_i = 2 * connec[i];
                        int dof_j = 2 * connec[j];
                        MatSetValues(A, 1, &dof_i,1, &dof_j, &matrix(2*i  ,2*j  ), ADD_VALUES);
                        
                        dof_i = 2 * connec[i] + 1;
                        dof_j = 2 * connec[j];
                        MatSetValues(A, 1, &dof_i, 1, &dof_j, &matrix(2*i+1,2*j  ), ADD_VALUES);
                        
                        dof_i = 2 * connec[i];
                        dof_j = 2 * connec[j] + 1;
                        MatSetValues(A, 1, &dof_i, 1, &dof_j, &matrix(2*i  ,2*j+1), ADD_VALUES);

                        dof_i = 2 * connec[i] + 1;
                        dof_j = 2 * connec[j] + 1;
                        MatSetValues(A, 1, &dof_i, 1, &dof_j, &matrix(2*i+1,2*j+1), ADD_VALUES);
                        
                        //Matrix Q and Qt
                        dof_i = 2 * connec[i];
                        dof_j = 2 * NNodes() + connec[j];
                        MatSetValues(A, 1, &dof_i, 1, &dof_j, &matrix(2*i  ,12+j), ADD_VALUES);
                    
                        dof_i = 2 * connec[i];
                        dof_j = 2 * NNodes() + connec[j];
                        MatSetValues(A, 1, &dof_j, 1, &dof_i, &matrix(12+j,2*i  ), ADD_VALUES);
                        
                        dof_i = 2 * connec[i] + 1;
                        dof_j = 2 * NNodes() + connec[j];
                        MatSetValues(A, 1, &dof_i, 1, &dof_j, &matrix(2*i+1,12+j), ADD_VALUES);
                        
                        dof_i = 2 * connec[i] + 1;
                        dof_j = 2 * NNodes() + connec[j];
                        MatSetValues(A, 1, &dof_j, 1, &dof_i, &matrix(12+j,2*i+1), ADD_VALUES);
                        
                        dof_i = 2 * NNodes() + connec[i];
                        dof_j = 2 * NNodes() + connec[j];
                        MatSetValues(A, 1, &dof_i, 1, &dof_j, &matrix(12+i,12+j), ADD_VALUES);
                    };
                    
                    //Rhs vector
                    int dof_i = 2 * connec[i];
                    VecSetValues(b, 1, &dof_i, &rhs[2*i  ], ADD_VALUES);
                    
                    dof_i = 2 * connec[i]+1;
                    VecSetValues(b, 1, &dof_i, &rhs[2*i+1], ADD_VALUES);
                    
                    dof_i = 2 * NNodes() + connec[i];
                    VecSetValues(b, 1, &dof_i, &rhs[12+i], ADD_VALUES);
                };
            };
        }; //Elements
        
        //Assemble matrices and vectors
        ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
        ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
        
        ierr = VecAssemblyBegin(b);
        ierr = VecAssemblyEnd(b);
        
        // MatView(A,PETSC_VIEWER_STDOUT_WORLD);
        // ierr = VecView(b,PETSC_VIEWER_STDOUT_WORLD);
        
        //Create KSP context to solve the linear system
        ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);
        
        ierr = KSPSetOperators(ksp,A,A);
        



    //     ierr = KSPSetTolerances(ksp,1.e-10,PETSC_DEFAULT,PETSC_DEFAULT,
    //                             500);
        
    //     ierr = KSPSetFromOptions(ksp);
        
    //     ierr = KSPGetPC(ksp,&pc);
        
    //     ierr = PCSetType(pc,PCJACOBI);
        
    //     //ierr = KSPSetType(ksp,KSPBCGS); 

    //     // ierr = KSPGMRESSetRestart(ksp, 10); 
        
    //        //ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
        

    //   //   ierr = MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,NULL,&nullsp);
    // // ierr = MatSetNullSpace(A, nullsp);
    // // ierr = MatNullSpaceDestroy(&nullsp);



#if defined(PETSC_HAVE_MUMPS)
        ierr = KSPSetType(ksp,KSPPREONLY);
        ierr = KSPGetPC(ksp,&pc);
        ierr = PCSetType(pc, PCLU);
#endif          
        ierr = KSPSetFromOptions(ksp);
        ierr = KSPSetUp(ksp);



        ierr = KSPSolve(ksp,b,u);

        ierr = KSPGetTotalIterations(ksp, &iterations);            

        //ierr = VecView(u,PETSC_VIEWER_STDOUT_WORLD);
        
        //Gathers the solution vector to the master process
        ierr = VecScatterCreateToAll(u, &ctx, &All);
        ierr = VecScatterBegin(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);
        ierr = VecScatterEnd(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);
        ierr = VecScatterDestroy(&ctx);

        ierr = VecScatterCreateToAll(b, &ctx, &Allu);
        ierr = VecScatterBegin(ctx, b, Allu, INSERT_VALUES, SCATTER_FORWARD);
        ierr = VecScatterEnd(ctx, b, Allu, INSERT_VALUES, SCATTER_FORWARD);
        ierr = VecScatterDestroy(&ctx);
        
        //Updates nodal values
        double p_;
        duNorm = 0.;
        double dpNorm = 0.;
        Ione = 1;
        
        for (int i = 0; i < NNodes(); ++i){
            Ii = 2*i;
            ierr = VecGetValues(All, Ione, &Ii, &val);
            NodeVec()[i] -> incrementAcceleration(0,val);
            NodeVec()[i] -> incrementVelocity(0,val*gamma*dTime);

            ierr = VecGetValues(Allu, Ione, &Ii, &val);
            duNorm += val*val;
        
            Ii = 2*i+1;
            ierr = VecGetValues(All, Ione, &Ii, &val);
            NodeVec()[i] -> incrementAcceleration(1,val);
            NodeVec()[i] -> incrementVelocity(1,val*gamma*dTime);

            ierr = VecGetValues(Allu, Ione, &Ii, &val);
            duNorm += val*val;
        };
        
        for (int i = 0; i<NNodes(); i++){
            Ii = 2*NNodes()+i;
            ierr = VecGetValues(All,Ione,&Ii,&val);
            p_ = val;
            NodeVec()[i] -> incrementPressure(p_);

            ierr = VecGetValues(Allu,Ione,&Ii,&val);
            dpNorm += val*val;
        };
        
        //Computes the solution vector norm
        //ierr = VecNorm(u,NORM_2,&val);

        std::clock_t t2 = std::clock();
     
        if(rank == 0){

            std::cout << "Iteration = " << inewton 
                      << " (" << iterations << ")"  
                      << "   Du Norm = " << std::scientific << sqrt(duNorm) 
                      << " " << sqrt(dpNorm)
                      << "  Time (s) = " << std::fixed
                      << 1000.*(t2-t1)/CLOCKS_PER_SEC/1000. << std::endl;
        };
                  
        ierr = KSPDestroy(&ksp); 
        ierr = VecDestroy(&b); 
        ierr = VecDestroy(&u); 
        ierr = VecDestroy(&All); 
        ierr = VecDestroy(&Allu); 
        ierr = MatDestroy(&A); 

        if (sqrt(duNorm) <= tolerance) {
            break;
        };
    };//Newton-Raphson
    
    return 0;
};
