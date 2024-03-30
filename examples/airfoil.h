
   // Defines the problem dimension
    const int dimension = 2;

{
    CompMesh* coarseModel = new CompMesh();
    TransientNavierStokes * matns = new TransientNavierStokes(17,2,1,0.1,0.01,0.9);
    // Stokes * matns = new Stokes(10,2,0.1,1.0);
    coarseModel->InsertMaterial(matns);
    
    //BC
    MatrixDouble val1(3,3);
    val1.setZero();
    VecDouble val2(3);
    VecDouble val3(3);
    val2.setZero();val3.setZero();


    val2[0] = 1.0;
    L2Projection * matbc3 = new L2Projection(18,2,0,val1,val2);
    val2.setZero();
    val2[1] = 1.;
    L2Projection * matbc1 = new L2Projection(19,2,3,val1,val2);
    
    L2Projection * matbc2 = new L2Projection(20,2,1,val1,val3);
    L2Projection * matbc4 = new L2Projection(21,2,0,val1,val3);

    coarseModel->InsertMaterial(matbc1);
    coarseModel->InsertMaterial(matbc2);
    coarseModel->InsertMaterial(matbc3);
    coarseModel->InsertMaterial(matbc4);

    GmshTools::Read(*coarseModel,"../airfoil.msh");

    TransientAnalysis an(coarseModel,SolverType::EUmfpack);
       
    std::vector<std::string> ScalarNames, VectorNames;
    ScalarNames = {"Pressure"};
    VectorNames = {"Velocity"};

    // an.PrintVariables("Dynamic",ScalarNames,VectorNames);
    an.Run();

    VTUGenerator::PrintResults(coarseModel,"resultCoarse",ScalarNames,VectorNames);

    
}           