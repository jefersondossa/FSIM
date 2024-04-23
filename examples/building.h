auto forcingFunctionPoisson = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];
    const auto &z=coord[2];
    force[0] = 1.-x*x;
    // force[0] = -2. * (x*x*(1.-3.*y) - (y-1.)*y*y + 3.*x*(y-1.)*y*y + x*x*x*(3.*y-1.));
};


auto exactSolPoisson = [](const VecDouble &coord, VecDouble &u, MatrixDouble &gradU){
    const auto &x=coord[0];
    const auto &y=coord[1];
    const auto &z=coord[2];

    // u[0] = x * x * (x-1.) * y * y * (y-1.);
    // gradU(0,0) = x*(3*x-2.)*(y-1.)*y*y;
    // gradU(1,0) = x*x*(x-1.)*y*(3.*y-2.);
    u[0] = z/1000.;
    // gradU(0,0) = 1.;
    // gradU(1,0) = 0.;
};


   // Defines the problem dimension
    const int dimension = 3;

{
    CompMesh* coarseModel = new CompMesh();
    // NavierStokes * matns = new NavierStokes(683,2,1,1.0);
    Stokes * matns = new Stokes(683,3,0.1,1.0);
    // Poisson * matns = new Poisson(683,3);
    // Elasticity3D * matns = new Elasticity3D(683,1.e3,.0);
    coarseModel->InsertMaterial(matns);
    
    //BC
    MatrixDouble val1(4,4);
    val1.setZero();
    VecDouble val2(4);
    VecDouble val3(4);
    val2.setZero();val3.setZero();
    //inflow
    val2[0] = 1.0;
    L2Projection * matbc1 = new L2Projection(677,3,0,val1,val2);
    matbc1->SetExactSolution(exactSolPoisson);
    val2.setZero();
    //outflow
    L2Projection * matbc2 = new L2Projection(678,3,1,val1,val2);

    //bottom and building
    L2Projection * matbc3 = new L2Projection(681,3,0,val1,val2);
    L2Projection * matbc4 = new L2Projection(682,3,0,val1,val2);
    
    val3[1] = 1.0;
    L2Projection * matbc5 = new L2Projection(679,2,3,val1,val3);
    val3.setZero();
    val3[2] = 1.;
    L2Projection * matbc6 = new L2Projection(680,2,3,val1,val3);
    // val3[0] = 1.0;
    // val3[1] = 1.0;
    // L2Projection * matbc4 = new L2Projection(9,2,3,val1,val3);
    coarseModel->InsertMaterial(matbc1);
    coarseModel->InsertMaterial(matbc2);
    coarseModel->InsertMaterial(matbc3);
    coarseModel->InsertMaterial(matbc4);
    coarseModel->InsertMaterial(matbc5);
    coarseModel->InsertMaterial(matbc6);

    GmshTools::Read(*coarseModel,"../TERRENO MALHA.msh");

    NonLinearAnalysis an(coarseModel,SolverType::ELDLt);

#ifdef RELEASE_BUILD
    // an.SType() = SolverType::EIterative;
    an.SType() = SolverType::EMumps;
    // an.SType() = SolverType::ECholmod;
#endif

    std::vector<std::string> ScalarNames, VectorNames;
    ScalarNames = {"Pressure"};
    VectorNames = {"Velocity"};

    // an.PrintVariables("Dynamic",ScalarNames,VectorNames);
    an.Run();

    VTUGenerator::PrintResults(coarseModel,"resultCoarse",ScalarNames,VectorNames);

    
}           