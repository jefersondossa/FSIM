//Analytic Solutions

auto exactSolPoisson = [](const VecDouble &coord, VecDouble &u, MatrixDouble &gradU){
    const auto &x=coord[0];
    const auto &y=coord[1];

    u[0] = x * x * (x-1.) * y * y * (y-1.);
    gradU(0,0) = x*(3*x-2.)*(y-1.)*y*y;
    gradU(1,0) = x*x*(x-1.)*y*(3.*y-2.);
    // u[0] = x * (x-1.) * y * (y-1.);
    // gradU(0,0) = (2*x-1.)*(y-1.)*y;
    // gradU(1,0) = x*(x-1.)*(2.*y-1.);
};

auto forcingFunctionPoisson = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];
    force[0] = -2. * (x*x*(1.-3.*y) - (y-1.)*y*y + 3.*x*(y-1.)*y*y + x*x*x*(3.*y-1.));
    // force[0] = -2.*(x-1.)*x - 2.*(y-1.)*y;
};


{

    GeoMesh * gmesh = new GeoMesh();
    GmshTools::Read(*gmesh,"../hierarquic2d.msh");

    gmesh->Print("gmesh.txt");

    CompMesh* cmesh = new CompMesh(gmesh, ApproxType::EHierarquic); 
    cmesh->SetDefaultOrder(2);
    // CompMesh* cmesh = new CompMesh(gmesh, ApproxType::EIsoparametric); 

    Poisson * matpoisson = new Poisson(8,2);
    cmesh->InsertMaterial(matpoisson);
    matpoisson->SetForcingFunction(forcingFunctionPoisson);
    matpoisson->SetExactSolution(exactSolPoisson);
    int nstate = 1;

    // Elasticity2D * matpoisson = new Elasticity2D(8,1.,.0);
    // coarseModel->InsertMaterial(matpoisson);
    // ElasticityPositional2D * matpoisson = new ElasticityPositional2D(8,1000.,.3);
    // coarseModel->InsertMaterial(matpoisson);
    // matpoisson->SetForcingFunction(forcingFunctionElasticity2D);
    // matpoisson->SetExactSolution(exactSolElasticity2D);
    //BC
    MatrixDouble val1(nstate,nstate);
    val1.setZero();
    VecDouble val2(nstate);
    VecDouble val3(nstate);
    val2.setZero();val3.setZero();
    // val2[0] = 1.0;
    L2Projection * matbc3 = new L2Projection(6,1,BoundaryConditionType::kDirichlet,val1,val2);
    val2.setZero();
    L2Projection * matbc1 = new L2Projection(5,1,BoundaryConditionType::kDirichlet,val1,val2);
    L2Projection * matbc2 = new L2Projection(7,1,BoundaryConditionType::kDirichlet,val1,val3);
    // matbc1->SetForcingFunction(forcingFunctionElasticity2D);
    // matbc1->SetExactSolution(exactSolElasticity2D);
    // matbc2->SetForcingFunction(forcingFunctionElasticity2D);
    // matbc2->SetExactSolution(exactSolElasticity2D);
    // matbc3->SetForcingFunction(forcingFunctionElasticity2D);
    // matbc3->SetExactSolution(exactSolElasticity2D);

    cmesh->InsertMaterial(matbc1);
    cmesh->InsertMaterial(matbc2);
    cmesh->InsertMaterial(matbc3);

    cmesh->AutoBuild();
    
    cmesh->Print("cmesh.txt");

    LinearAnalysis an(cmesh,SolverType::ELDLt);
    an.Run();
    // an.PrintGlobalMatrix();
    // an.PrintGlobalRhs();

    std::vector<std::string> ScalarNames, VectorNames;
    ScalarNames = {"Solution","ExactSolution"};
    // VectorNames = {"Displacement"};

    VTUGenerator::PrintResults(cmesh,"result2d",ScalarNames,VectorNames);

    // VecDouble errors;
    // an.PostProcessError(errors);

}           