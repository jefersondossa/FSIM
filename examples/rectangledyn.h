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
    // force[0] = 1.-x*x;
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
    const int degree = 1;

{
// //  Create problem variables 
    
    Arlequin   arlequinProblem; 

//==========================================================================
//==============================PROBLEM MESHES==============================
//==========================================================================
   
    // ProblemParameters pParameters;
    // // pParameters.ProbType() = ESolidPositional;
    // pParameters.ProbType() = EPoisson;
    // if (pParameters.ProbType() == EPoisson){
    //     pParameters.setForcingFunction(forcingFunctionPoisson);
    //     pParameters.setExactSolution(exactSolPoisson);
    // } else if (pParameters.ProbType() == EElastic || pParameters.ProbType() == ESolidPositional){
    //     pParameters.SetElasticity(1.,0.);
    //     // pParameters.setForcingFunction(forcingFunctionElasticity2D);
    //     // pParameters.setExactSolution(exactSolElasticity2D);
    // } else if (pParameters.ProbType() == EStokes){
    //     pParameters.SetIncompressibleFluid(0.01,1.);
    //     pParameters.setTimeStep(1.);//Needed for stabilization parameter
    //     pParameters.setForcingFunction(forcingFunctionStokes);
    //     pParameters.setExactSolution(exactSolStokes);
    // } else if (pParameters.ProbType() == ENavierStokes){
    //     pParameters.SetIncompressibleFluid(0.01,1.);
    //     pParameters.setTimeStep(1.);//Needed for stabilization parameter
    //     pParameters.setForcingFunction(forcingFunctionNavierStokes);
    //     pParameters.setExactSolution(exactSolStokes);
    // }
    // pParameters.Solver() = SolverType::EUmfpack;
    // pParameters.ArlequinStab() = ArlequinStabType::ENoStab;
    // // pParameters.ArlequinStab() = ArlequinStabType::EOption1;
    // pParameters.setSpectralRadius(1.);
    // pParameters.setArlequinOperatorConstants(1.,0.0);

    CompMesh* coarseModel = new CompMesh();
    // CompMesh* fineModel = new CompMesh();  

    // Poisson * matpoisson = new Poisson(8,2);
    // coarseModel->InsertMaterial(matpoisson);
    // matpoisson->SetForcingFunction(forcingFunctionPoisson);
    // matpoisson->SetExactSolution(exactSolPoisson);

    TransientElasticity2D * matpoisson = new TransientElasticity2D(8,1.,.0,true,0.,1.,0.1);
    coarseModel->InsertMaterial(matpoisson);
    
    //BC
    MatrixDouble val1(2,2);
    val1.setZero();
    VecDouble val2(2);
    VecDouble val3(2);
    val2.setZero();val3.setZero();
    val2[0] = 1.0;
    L2Projection * matbc3 = new L2Projection(6,1,3,val1,val2);
    val2.setZero();
    L2Projection * matbc1 = new L2Projection(5,1,0,val1,val2);
    L2Projection * matbc2 = new L2Projection(7,1,1,val1,val3);
    coarseModel->InsertMaterial(matbc1);
    coarseModel->InsertMaterial(matbc2);
    coarseModel->InsertMaterial(matbc3);

    GmshTools::Read(*coarseModel,"../coarse_test.msh");
    // GmshTools::Read(*fineModel,"../fine_test.msh");

    TransientAnalysis an(coarseModel,SolverType::ELDLt,true);
       
    std::vector<std::string> ScalarNames, VectorNames;
    ScalarNames = {"SigmaX","ExactSigmaX"};
    VectorNames = {"Displacement","ExactDisplacement"};

    an.PrintVariables("Dynamic",ScalarNames,VectorNames);
    an.Run(100);
    
}           