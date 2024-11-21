
   // Defines the problem dimension
    const int dimension = 2;

{
    CompMesh* coarseModel = new CompMesh();
    TransientNavierStokes * matns = new TransientNavierStokes(17,2,1.,1.e0,1.,0.0);
    // Stokes * matns = new Stokes(17,2,1.,100.0);
    // NavierStokes * matns = new NavierStokes(17,2,1.,1.0);

    coarseModel->InsertMaterial(matns);
    
    int nstate = matns->NState();
    //BC
    MatrixDouble val1(nstate,nstate);
    val1.setZero();
    VecDouble val2(nstate);
    VecDouble val3(nstate);
    val2.setZero();val3.setZero();


    //Inflow
    val2[0] = 1.0;
    val2[1] = 0.00001;
    L2Projection * matbc3 = new L2Projection(18,2,BoundaryConditionType::kDirectionalNonHomogeneousDirichlet,val1,val2);
    //Up and Down
    val2.setZero();
    val2[0] = 1.;
    val2[1] = 0.00001;
    L2Projection * matbc1 = new L2Projection(19,2,BoundaryConditionType::kDirectionalNonHomogeneousDirichlet,val1,val2);
    
    //Outflow
    L2Projection * matbc2 = new L2Projection(20,2,BoundaryConditionType::kNeumann,val1,val3);
    //Airfoil
    val3[0] = 1.0;
    val3[1] = 1.0;
    L2Projection * matbc4 = new L2Projection(21,2,BoundaryConditionType::kDirectionalHomogeneousDirichlet,val1,val3);

    coarseModel->InsertMaterial(matbc1);
    coarseModel->InsertMaterial(matbc2);
    coarseModel->InsertMaterial(matbc3);
    coarseModel->InsertMaterial(matbc4);

    GmshTools::Read(*coarseModel,"../airfoil.msh");

    VecDouble InitSol(nstate);
    InitSol.setZero();
    InitSol[0] = 1.0;
    coarseModel->SetSolution(InitSol);

    TransientAnalysis an(coarseModel,SolverType::EMumps,false);
    // NonLinearAnalysis an(coarseModel,SolverType::EMumps,1.e-6,10);
    // LinearAnalysis an(coarseModel,SolverType::EMumps);
       
    std::vector<std::string> ScalarNames, VectorNames;
    ScalarNames = {"Pressure"};
    VectorNames = {"Velocity","Acceleration"};
    // VectorNames = {"Displacement"};
    VTUGenerator::PrintResults(coarseModel,"resultCoarse",ScalarNames,VectorNames);

    an.PrintVariables("Airfoil",ScalarNames,VectorNames);
    an.Run(100);
    // an.Run();

    // VTUGenerator::PrintResults(coarseModel,"resultCoarse",ScalarNames,VectorNames);

    
}           