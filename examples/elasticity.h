// Defines the problem dimension
const int dimension = 2;
auto forcingFunction = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];
    force[0] = 10E3;
};

{
    CompMesh* cmesh = new CompMesh();
    Elasticity2D * matelasticity = new Elasticity2D(3, 2E11, 0.3);
    //matpoisson->SetForcingFunction(forcingFunctionPoisson);
    cmesh->InsertMaterial(matelasticity);
    
    //BC
    MatrixDouble val1(2,2);
    val1.setZero();
    VecDouble val2(2);
    val2.setZero();
    L2Projection * matbc1 = new L2Projection(2,dimension-1,BoundaryConditionType::kDirichlet,val1,val2);
    val2[0] = 10E3;
    L2Projection * matbc2 = new L2Projection(1,dimension-2,BoundaryConditionType::kNeumann,val1,val2);
    cmesh->InsertMaterial(matbc1);
    cmesh->InsertMaterial(matbc2);
    
    GmshTools::Read(*cmesh,"../teste.msh");

    LinearAnalysis an(cmesh,SolverType::ELDLt);
       
    std::vector<std::string> ScalarNames, VectorNames;
    ScalarNames = {"SigmaX","SigmaY","TauXY"};
    VectorNames = {"Displacement"};

    an.Run();

    VTUGenerator::PrintResults(cmesh,"result",ScalarNames,VectorNames);    
    //VTUGenerator::PrintResultsGraph(cmesh,"result",ScalarNames,VectorNames);    
}   