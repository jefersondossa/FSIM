   // Defines the problem dimension
    const int dimension = 2;
    const int degree = 1;

{
// //  Create problem variables 
    
    Arlequin   arlequinProblem; 

//==========================================================================
//==============================PROBLEM MESHES==============================
//==========================================================================
   

    CompMesh* coarseModel = new CompMesh();
    CompMesh* fineModel = new CompMesh();  

    //COARSE MODEL
    Elasticity2D * matelas = new Elasticity2D(14,1.e10,0.0,true);
    coarseModel->InsertMaterial(matelas);
    
    //BC
    //Load
    MatrixDouble val1(2,2);
    val1.setZero();
    VecDouble val2(2);
    val2.setZero();
    val2[1] = -10.;
    L2Projection * matbc1 = new L2Projection(11,2,1,val1,val2);
    
    //Left Support
    val2.setZero();
    L2Projection * matbc2 = new L2Projection(12,2,0,val1,val2);
    val2[1] = 1.;
    L2Projection * matbc3 = new L2Projection(13,2,3,val1,val2);
    
    coarseModel->InsertMaterial(matbc1);
    coarseModel->InsertMaterial(matbc2);
    coarseModel->InsertMaterial(matbc3);

    //FINE MODEL
    Elasticity2D * matelas2 = new Elasticity2D(16,1.e10,0.0,true);
    fineModel->InsertMaterial(matelas2);
    //BC
    //Load
    val2[1] = -10.; 
    L2Projection * matbc4 = new L2Projection(15,2,1,val1,val2);
    fineModel->InsertMaterial(matbc4);
    //Glue
    val2.setZero();
    L2Projection * matbc5 = new L2Projection(19,2,1,val1,val2);
    fineModel->InsertMaterial(matbc5);

    GmshTools::Read(*coarseModel,"../BeamCoarse.msh");
    GmshTools::Read(*fineModel,"../BeamFine.msh");

    std::set<int> gluematids={19};
    std::vector<CompMesh *> meshvector(2);
    meshvector[0] = coarseModel;
    meshvector[1] = fineModel;
    Arlequin arl(meshvector,1.,0.);
    arl.SetGlueZoneThichkess(0.45);
    arl.SetGlueIds(gluematids);
    std::vector<std::string> ScalarNames, VectorNames;
    VTUGenerator::PrintResults(coarseModel,"resultCoarse",ScalarNames,VectorNames);
    VTUGenerator::PrintResults(fineModel,"resultFine",ScalarNames,VectorNames);
    // VTUGenerator::PrintResults(&arl,"ArlequinMesh");

    arl.SetUp();

    

    // LinearAnalysis an(coarseModel,SolverType::EUmfpack);
    // an.Run();
    LinearAnalysis an(&arl,SolverType::ELDLt);
    // NonLinearAnalysis an(&arl,SolverType::EUmfpack,1.e-6,2);
    // LinearAnalysis an(coarseModel,SolverType::ELDLt);
    // NonLinearAnalysis an(coarseModel,SolverType::ELDLt);
    // NonLinearAnalysis an(coarseModel,SolverType::EKLU);
    // NonLinearAnalysis an(coarseModel,SolverType::ESPQR);
    an.Run();

    ScalarNames = {};
    VectorNames = {"Displacement","Stress","Strain"};

    VTUGenerator::PrintResults(coarseModel,"resultCoarse",ScalarNames,VectorNames);
    VTUGenerator::PrintResults(fineModel,"resultFine",ScalarNames,VectorNames);
    // VTUGenerator::PrintResults(arl.MeshVec()[2],"resultCoupling");
}           