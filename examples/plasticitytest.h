
   // Defines the problem dimension
    const int dimension = 2;
{

    CompMesh* cmesh = new CompMesh();

    Elasticity2D * matelas = new Elasticity2D(8,1000.,.0);
    //BC
    MatrixDouble val1(2,2);
    val1.setZero();
    VecDouble val2(2);
    VecDouble val3(2);
    val2.setZero();val3.setZero();
    val2[0] = 1.5;
    L2Projection * matbc3 = new L2Projection(6,1,0,val1,val2);
    val2.setZero();
    L2Projection * matbc1 = new L2Projection(5,1,0,val1,val2);
    L2Projection * matbc2 = new L2Projection(7,1,1,val1,val3);
    

    cmesh->InsertMaterial(matbc1);
    cmesh->InsertMaterial(matbc2);
    cmesh->InsertMaterial(matbc3);
    LinearHardening *plastmodel = new LinearHardening(matelas,111.,100.); 
    cmesh->InsertMaterial(plastmodel);
    std::vector<L2Projection *> bcIncrement = {matbc3};
  
    GmshTools::Read(*cmesh,"../coarse_test.msh");

    IncrementalAnalysis an(cmesh,SolverType::ESuiteSparse, 10, bcIncrement);
    
    std::vector<std::string> ScalarNames, VectorNames;
    ScalarNames = {"Stress"};
    VectorNames = {"Displacement"};
    an.Run("plasticitytest",ScalarNames,VectorNames);

    
}           