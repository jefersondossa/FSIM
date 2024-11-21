   // Defines the problem dimension
    const int dimension = 2;

auto forcingFunctionPoisson = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];
    force[0] = 2.*(1.-coord[0])*coord[0]+2.*(1-coord[1])*coord[1];
    // force[0] = 1.-x*x;
    // force[0] = -2. * (x*x*(1.-3.*y) - (y-1.)*y*y + 3.*x*(y-1.)*y*y + x*x*x*(3.*y-1.));
};

{
    CompMesh* cmesh = new CompMesh();
    Poisson * matpoisson = new Poisson(1,dimension);
    matpoisson->SetForcingFunction(forcingFunctionPoisson);
    cmesh->InsertMaterial(matpoisson);
    
    //BC
    MatrixDouble val1(1,1);
    val1.setZero();
    VecDouble val2(1);
    val2.setZero();
    L2Projection * matbc1 = new L2Projection(2,dimension-1,0,val1,val2);
    val2[0] = 0.0;
    L2Projection * matbc2 = new L2Projection(3,dimension-2,1,val1,val2);
    cmesh->InsertMaterial(matbc1);
    cmesh->InsertMaterial(matbc2);
    
    GmshTools::Read(*cmesh,"../quads.msh");

    LinearAnalysis an(cmesh,SolverType::ELU);
       
    std::vector<std::string> ScalarNames, VectorNames;
    ScalarNames = {"Solution"};
    VectorNames = {"Derivative"};

    an.Run();

    // VTUGenerator::PrintResults(cmesh,"result",ScalarNames,VectorNames);    
    VTUGenerator::PrintResultsGraph(cmesh,"result",ScalarNames,VectorNames);    
}           