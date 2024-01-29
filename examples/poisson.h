   // Defines the problem dimension
    const int dimension = 1;

{
    CompMesh* cmesh = new CompMesh();
    TransientPoisson * matpoisson = new TransientPoisson(4,dimension,0.,1.5,0.05);
    cmesh->InsertMaterial(matpoisson);
    
    //BC
    MatrixDouble val1(1,1);
    val1.setZero();
    VecDouble val2(1);
    val2.setZero();
    L2Projection * matbc1 = new L2Projection(2,dimension-1,0,val1,val2);
    val2[0] = 1.0;
    L2Projection * matbc2 = new L2Projection(3,dimension-1,1,val1,val2);
    cmesh->InsertMaterial(matbc1);
    cmesh->InsertMaterial(matbc2);
    
    GmshTools::Read(*cmesh,"../line.msh");

    TransientAnalysis an(cmesh,SolverType::ELDLt,true);
       
    std::vector<std::string> ScalarNames, VectorNames;
    ScalarNames = {"Solution"};
    VectorNames = {"Derivative"};

    an.PrintVariables("Dynamic",ScalarNames,VectorNames);
    an.Run(100);
    
}           