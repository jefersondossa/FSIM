//Analytic Solutions

auto exactSolPoisson = [](const VecDouble &coord, VecDouble &u, MatrixDouble &gradU){
    const auto &x=coord[0];
    const auto &y=coord[1];

    // u[0] = x * x * (x-1.) * y * y * (y-1.);
    // gradU(0,0) = x*(3*x-2.)*(y-1.)*y*y;
    // gradU(1,0) = x*x*(x-1.)*y*(3.*y-2.);
    u[0] = x;
    gradU(0,0) = 1.;
    gradU(1,0) = 0.;
};

auto forcingFunctionPoisson = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];
    force[0] = 1.-x*x;
    // force[0] = -2. * (x*x*(1.-3.*y) - (y-1.)*y*y + 3.*x*(y-1.)*y*y + x*x*x*(3.*y-1.));
};

auto exactSolElasticity2D = [](const VecDouble &coord, VecDouble &u, MatrixDouble &gradU){
    const auto &x=coord[0];
    const auto &y=coord[1];
    auto pi = M_PI;
    u[0] = cos(pi*x)*sin(2.*pi*y);
    u[1] = cos(pi*y)*sin(pi*x);
    gradU(0,0) = -pi*sin(pi*x)*sin(2.*pi*y);
    gradU(0,1) = pi*cos(pi*x)*cos(pi*y);
    gradU(1,0) = 2.*pi*cos(pi*x)*cos(2.*pi*y);
    gradU(1,1) = -pi*sin(pi*x)*sin(pi*y);
    // u[0] = x;
    // u[1] = 0.;
    // gradU(0,0) = 1.;
    // gradU(0,1) = 0.;
    // gradU(1,0) = 0.;
    // gradU(1,1) = 0.;
};

auto forcingFunctionElasticity2D = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];
    double E=1.;
    double poisson=0.;
    auto pi = M_PI;
    force[0] = -(E*pi*pi*cos(pi*x)*((1.+poisson)*sin(pi*y) + 2.*(3. - 2.*poisson)*sin(2.*pi*y)))/(2.*(-1. + poisson*poisson));
    force[1] = -(E*pi*pi*(-((-3. + poisson)*cos(pi*y)) + 2.*(1. + poisson)*cos(2*pi*y))*sin(pi*x))/(2.*(-1. + poisson*poisson));
};


{
    CompMesh* cmesh = new CompMesh(ApproxType::EHierarquic); 
    cmesh->SetDefaultOrder(2);
    // CompMesh* cmesh = new CompMesh(ApproxType::EIsoparametric); 

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

    GmshTools::Read(*cmesh,"../hierarquic2d.msh");
    cmesh->Print("cmesh.txt");

    LinearAnalysis an(cmesh,SolverType::ELDLt);
    an.Run();

    std::vector<std::string> ScalarNames, VectorNames;
    ScalarNames = {"Solution"};
    // VectorNames = {"Displacement"};

    VTUGenerator::PrintResults(cmesh,"result2d",ScalarNames,VectorNames);

    // VecDouble errors;
    // an.PostProcessError(errors);

}           