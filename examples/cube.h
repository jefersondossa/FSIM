//Analytic Solutions
   // Defines the problem dimension
    const int dimension = 3;
{
// //  Create problem variables 
    
//==========================================================================
//==============================PROBLEM MESHES==============================
//==========================================================================
auto yieldFunction = [](const double &plast, double &sigma_y, double &hardening){
    // sigma_y =.5/sqrt(3.)+500.*plast;
    // sigma_y =50.+500.*plast;
    // hardening = 500.;

    // sigma_y = 20e3+5000000.000*plast;
    // // sigma_y = 10.0;//+0.0001*plast;
    // hardening = 5000000.000;
    sigma_y = .5+500.00*plast;

    hardening = 500.00;
};

    CompMesh* coarseModel = new CompMesh();

    // Elasticity3D * matpoisson = new Elasticity3D(15,1.e3,.0);
    Poisson * matpoisson = new Poisson(15,3,1);
    coarseModel->InsertMaterial(matpoisson);
    // VonMises *plastmodel = new VonMises(matpoisson);
    // DruckerPrager *plastmodel = new DruckerPrager(matpoisson,00.,0.,false);
    // plastmodel->SetUniaxialYieldFunction(yieldFunction);
    // coarseModel->InsertMaterial(plastmodel);

    //BC
    MatrixDouble val1(1,1);
    val1.setZero();
    VecDouble val2(1);
    val2.setZero();
    val2[0] = 0.;
    L2Projection * matbc1 = new L2Projection(16,3,0,val1,val2);
    val2.setZero();
    // L2Projection * matbc3 = new L2Projection(18,3,1,val1,val2);
    // val2[0] = -20.0;
    val2[0] = 1.;
    // val2[2] = -20.0;
    // val2[1] = 0.4;
    L2Projection * matbc2 = new L2Projection(17,3,0,val1,val2);
    val2.setZero();
    // val2[0] = 1.;
    L2Projection * matbc4 = new L2Projection(18,3,1,val1,val2);
    val2.setZero();
    // val2[0] = 1.;
    L2Projection * matbc5 = new L2Projection(19,3,1,val1,val2);
    // val2.setZero();
    // val2[2] = 1.;
    // L2Projection * matbc6 = new L2Projection(21,3,3,val1,val2);
    // matbc3->SetForcingFunction(forcingFunctionElasticity3D);
    // matbc3->SetExactSolution(exactSolElasticity3D);

    coarseModel->InsertMaterial(matbc1);
    coarseModel->InsertMaterial(matbc2);
    // coarseModel->InsertMaterial(matbc3);
    coarseModel->InsertMaterial(matbc4);
    coarseModel->InsertMaterial(matbc5);
    // coarseModel->InsertMaterial(matbc6);

    GmshTools::Read(*coarseModel,"../cube.msh");

    // NonLinearAnalysis an(coarseModel,SolverType::ELU);
    // LinearAnalysis an(coarseModel,SolverType::ELU);

    std::vector<L2Projection *> bcIncrement = {matbc2};
    NonLinearAnalysis an(coarseModel,SolverType::ELU, 1.e-7,100);
#ifdef RELEASE_BUILD
    an.SType() = SolverType::EUmfpack;
    // an.SType() = SolverType::ECholmod;
#endif    
    VecDouble increment(3);
    // increment.setZero();
    increment[0] = 0.000035;
    // an.SetIncrement(increment);
    std::vector<std::string> ScalarNames, VectorNames;
    // ScalarNames = {"PlasticStrain"};
    ScalarNames = {"Solution"};
    // VectorNames = {"Derivative"};
    // VectorNames = {"Displacement","Stress"};
    an.Run();
    // an.Run("pureShear",ScalarNames,VectorNames);

    // std::vector<std::string> ScalarNames, VectorNames;
    // ScalarNames = {};
    // VectorNames = {"Displacement"};
    // // VectorNames = {"Displacement","ExactDisplacement"};
    // an.GlobalMatrix()->ZeroSolution();
    // int64_t row = 9;
    // double val = 1.;
    // an.GlobalMatrix()->PutValueSolution(row,val);
    // an.UpdateSolution();

    VTUGenerator::PrintResults(coarseModel,"cube",ScalarNames,VectorNames);
    // VTUGenerator::PrintResults(fineModel,"resultFine");
    // VTUGenerator::PrintResults(arl.MeshVec()[2],"resultCoupling");

    // VecDouble errors;
    // an.PostProcessError(errors);

}           