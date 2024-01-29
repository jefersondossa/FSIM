//Analytic Solutions

auto exactSolElasticity3D = [](const VecDouble &coord, VecDouble &u, MatrixDouble &gradU){
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

auto forcingFunctionElasticity3D = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];
    double E=1.;
    double poisson=0.;
    auto pi = M_PI;
    force[0] = -(E*pi*pi*cos(pi*x)*((1.+poisson)*sin(pi*y) + 2.*(3. - 2.*poisson)*sin(2.*pi*y)))/(2.*(-1. + poisson*poisson));
    force[1] = -(E*pi*pi*(-((-3. + poisson)*cos(pi*y)) + 2.*(1. + poisson)*cos(2*pi*y))*sin(pi*x))/(2.*(-1. + poisson*poisson));
};


   // Defines the problem dimension
    const int dimension = 3;
{
// //  Create problem variables 
    
//==========================================================================
//==============================PROBLEM MESHES==============================
//==========================================================================
   
	MPI_Barrier(PETSC_COMM_WORLD);   

    CompMesh* coarseModel = new CompMesh();

    Elasticity3D * matpoisson = new Elasticity3D(15,1.,.0);
    coarseModel->InsertMaterial(matpoisson);
    // matpoisson->SetForcingFunction(forcingFunctionElasticity3D);
    // matpoisson->SetExactSolution(exactSolElasticity3D);
    //BC
    MatrixDouble val1(3,3);
    val1.setZero();
    VecDouble val2(3);
    VecDouble val3(3);
    val2.setZero();val3.setZero();
    val2[0]=1.;
    L2Projection * matbc3 = new L2Projection(13,3,0,val1,val2);
    // matbc3->SetForcingFunction(forcingFunctionElasticity3D);
    // matbc3->SetExactSolution(exactSolElasticity3D);

    coarseModel->InsertMaterial(matbc3);

    GmshTools::Read(*coarseModel,"../cube.msh");

    // NonLinearAnalysis an(coarseModel,SolverType::ECholmod);
    LinearAnalysis an(coarseModel,SolverType::ECholmod);
    an.Run();

    std::vector<std::string> ScalarNames, VectorNames;
    ScalarNames = {};
    VectorNames = {"Displacement"};
    // VectorNames = {"Displacement","ExactDisplacement"};

    VTUGenerator::PrintResults(coarseModel,"resultCoarse",ScalarNames,VectorNames);
    // VTUGenerator::PrintResults(fineModel,"resultFine");
    // VTUGenerator::PrintResults(arl.MeshVec()[2],"resultCoupling");

    // VecDouble errors;
    // an.PostProcessError(errors);

}           