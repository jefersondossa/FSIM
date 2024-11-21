auto exactSolStokes = [](const VecDouble &coord, VecDouble &u, MatrixDouble &gradU){
    const auto &x=coord[0];
    const auto &y=coord[1];

    u[0] = -sin(x)*sin(y);
    u[1] = -cos(x)*cos(y);
    u[2] = cos(x)*cos(y);//-sin(2.*pi*x)*cos(2.*pi*y);
    gradU(0,0) = -cos(x)*sin(y);
    gradU(0,1) = -sin(x)*cos(y);
    gradU(1,0) =  sin(x)*cos(y);
    gradU(1,1) =  cos(x)*sin(y);
    gradU(0,2) = -sin(x)*sin(y);
    gradU(1,2) =  cos(x)*cos(y);
};

auto forcingFunctionStokes = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];
    double visc = 0.01;

    force[0] = -(1+2.*visc)*sin(x)*sin(y);
    force[1] = (1-2.*visc)*cos(x)*cos(y);
};

auto forcingFunctionNavierStokes = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];
    double visc = 0.01;

    force[0] = -(sin(x)*(cos(x - 2.*y) + cos(x + 2.*y) - 2.*(1. + 2.*visc)*sin(y)))/2.;
    force[1] = -cos(y)*(cos(x) - 2.*visc*cos(x) + cos(x)*cos(x)*sin(y) - sin(x)*sin(x)*sin(y));
};
   // Defines the problem dimension
    const int dimension = 2;

{
    CompMesh* coarseModel = new CompMesh();
    NavierStokes * matns = new NavierStokes(10,2,1,1.e1);
    // Poisson * matns = new Poisson(10,2,1);
    // Stokes * matns = new Stokes(10,2,0.1,1.0);
    coarseModel->InsertMaterial(matns);
    
    int nstate = matns->NState();
    //BC
    MatrixDouble val1(nstate,nstate);
    val1.setZero();
    VecDouble val2(nstate);
    VecDouble val3(nstate);
    val2.setZero();val3.setZero();
    val2[0] = 1.0;
    val2[1] = 0.000000001;
    L2Projection * matbc3 = new L2Projection(6,2,BoundaryConditionType::kDirectionalNonHomogeneousDirichlet,val1,val2);
    val2.setZero();
    // val2[2]=1.;
    L2Projection * matbc1 = new L2Projection(7,2,BoundaryConditionType::kNeumann,val1,val2);
    val3[1] = 1.0;
    val3[0] = 1.0;
    L2Projection * matbc2 = new L2Projection(8,2,BoundaryConditionType::kDirectionalHomogeneousDirichlet,val1,val3);
    // val3[0] = 1.0;
    // val3[1] = 1.0;
    // L2Projection * matbc4 = new L2Projection(9,2,3,val1,val3);
    coarseModel->InsertMaterial(matbc1);
    coarseModel->InsertMaterial(matbc2);
    coarseModel->InsertMaterial(matbc3);
    // coarseModel->InsertMaterial(matbc4);

    GmshTools::Read(*coarseModel,"../cylinder.msh");

    NonLinearAnalysis an(coarseModel,SolverType::EUmfpack,1.e-6,1);
       
    std::vector<std::string> ScalarNames, VectorNames;
    ScalarNames = {"Pressure"};
    VectorNames = {"Velocity"};

    // an.PrintVariables("Dynamic",ScalarNames,VectorNames);
    an.Run();

    VTUGenerator::PrintResults(coarseModel,"resultCoarse",ScalarNames,VectorNames);

    
}           