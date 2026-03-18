// Defines the problem dimension
const int dimension = 2;
auto forcingFunction = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];
    force[0] = -1*y;
};

{
    //GLOBAL MODEL
    CompMesh* cmeshG = new CompMesh();
    Elasticity2D* matelasticityG = new Elasticity2D(1, 1.0, 0.0);
    cmeshG->InsertMaterial(matelasticityG);
    
    //BC 
    MatrixDouble val1(2,2);
    val1.setZero();
    VecDouble val2(2);
    val2.setZero();
    L2Projection * matbcG1 = new L2Projection(2,dimension-1,BoundaryConditionType::kDirichlet,val1,val2);
    val2[1] = 1.0;
    L2Projection * matbcG2 = new L2Projection(3,dimension-1,BoundaryConditionType::kDirectionalHomogeneousDirichlet,val1,val2);
    L2Projection * matbcG3 = new L2Projection(4,dimension-1,BoundaryConditionType::kNeumann,val1,val2);
    matbcG3->SetForcingFunction(forcingFunction);

    cmeshG->InsertMaterial(matbcG1);
    cmeshG->InsertMaterial(matbcG2);
    cmeshG->InsertMaterial(matbcG3);
    
    GmshTools::Read(*cmeshG,"../global.msh");

    LinearAnalysis anG(cmeshG,SolverType::ELDLt);
       
    std::vector<std::string> ScalarNames, VectorNames;
    ScalarNames = {"SigmaX","SigmaY","TauXY"};
    VectorNames = {"Displacement"};

    anG.Run();

    VTUGenerator::PrintResults(cmeshG,"globalResult",ScalarNames,VectorNames);    
    //VTUGenerator::PrintResultsGraph(cmesh,"result",ScalarNames,VectorNames);    

    //LOCAL MODEL
    CompMesh* cmeshL = new CompMesh();
    Elasticity2D* matelasticityL = new Elasticity2D(1, 1.0, 0.0);
    cmeshL->InsertMaterial(matelasticityL);

    //BC
    val1.setZero();
    val2.setZero();
    val2[1] = 1.0;
    L2Projection * matbcL1 = new L2Projection(3,dimension-1,BoundaryConditionType::kDirectionalHomogeneousDirichlet,val1,val2);
    val2.setZero();
    val2[0] = -2.0;
    L2Projection * matbcL2 = new L2Projection(2,dimension-1,BoundaryConditionType::kDirectionalNonHomogeneousDirichlet,val1,val2);
    L2Projection * matbcL3 = new L2Projection(4,dimension-1,BoundaryConditionType::kNeumann,val1,val2);
    matbcL3->SetForcingFunction(forcingFunction);

    cmeshL->InsertMaterial(matbcL1);
    cmeshL->InsertMaterial(matbcL2);
    cmeshL->InsertMaterial(matbcL3);
    
    GmshTools::Read(*cmeshL,"../local.msh");

    LinearAnalysis anL(cmeshL,SolverType::ELDLt);

    anL.Run();
    //anL.PrintGlobalMatrix();

    VTUGenerator::PrintResults(cmeshL,"localResult",ScalarNames,VectorNames); 
}   