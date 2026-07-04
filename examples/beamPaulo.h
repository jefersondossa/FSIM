
   // Defines the problem dimension
    const int dimension = 2;
{
    CompMesh* cmesh = new CompMesh(); 

    // Elasticity2D * matelas = new Elasticity2D(11,210.e9,.3,true,0.05);
    // ElasticityPositional2D * matelas = new ElasticityPositional2D(13,10.e6,.0,0.1,true);
    Elasticity2D * matelas = new Elasticity2D(13,10.e6,.0,0.1,true);

    //BC
    MatrixDouble val1(2,2);
    val1.setZero();
    VecDouble val2(2);
    val2.setZero();
    //Encastre
    // val2[0] = 1.;
    L2Projection * matbc1 = new L2Projection(11,2,BoundaryConditionType::kDirichlet,val1,val2);

    val2.setZero();
    //Encastre point
    val2[1] = -5.4e3;
    L2Projection * matbc2 = new L2Projection(12,2,BoundaryConditionType::kNeumann,val1,val2);

    // //Load
    // val2[1] = 000.;
    // L2Projection * matbc3 = new L2Projection(9,2,1,val1,val2);
    
    cmesh->InsertMaterial(matbc1);
    cmesh->InsertMaterial(matbc2);
    // cmesh->InsertMaterial(matbc3);
    
    // VonMises *plastmodel = new VonMises(matelas);
    // plastmodel->SetUniaxialYieldFunction(yieldFunction);
    cmesh->InsertMaterial(matelas);
    // cmesh->InsertMaterial(plastmodel); 
    std::vector<L2Projection *> bcIncrement = {matbc2};
   
    GmshTools::Read(*cmesh,"../dam.msh");

    IncrementalAnalysis an(cmesh,SolverType::EUmfpack, 10, bcIncrement,1.e-7,300);

    // NonLinearAnalysis an(cmesh,SolverType::ELU, 1.e-7,300);
    // VecDouble increment(2);
    // increment.setZero();    
    // increment[1] = -1000.;
    // increment[1] = 0.02;
    // an.SetIncrement(increment);
#ifdef RELEASE_BUILD
    an.SType() = SolverType::EUmfpack;
    // an.SType() = SolverType::ECholmod;
#endif

    
    // VectorNames = {"Displacement","Stress","Strain"};
    // an.Run();
    
    std::vector<std::string> ScalarNames, VectorNames;
    VectorNames = {"Displacement","Stress"};
    // VTUGenerator::PrintResults(cmesh,"testPaulo",ScalarNames,VectorNames);
    an.Run("testPaulo",ScalarNames,VectorNames);

    
    // std::vector<std::string> integrate = {"Solution","DerivativeX","DerivativeY"};
    // std::set<int> matid = {6};
    // std::map<std::string,VecDouble> result;
    // cmesh->Integrate(matid,integrate,result);
    // REAL aux = 206.9/((1.29)*(1.-2.*.29));
    // REAL valX = aux * ((1.-.29)*result["DerivativeX"][0] + .29*result["DerivativeY"][0])/.45;
    // REAL valY = aux * (.29*result["DerivativeX"][0] + (1.-.29)*result["DerivativeY"][0])/.45;
    // int a = 0;
    
}           