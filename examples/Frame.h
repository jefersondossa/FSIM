
   // Defines the problem dimension
const int dimension = 1;
auto forcing = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];
    force[1] = x-1.;
};
auto forcingTransient = [](const VecDouble &coord, double &t, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];
    force[1] = 100.;//*std::sin(0.2*t);
};


{

    CompMesh* cmesh = new CompMesh(); 

    // LinearFrame * bhorizontal = new LinearFrame(5,1.,1./12.,1.);
    TransientPositionalFrame2D * bhorizontal = new TransientPositionalFrame2D(5,25.e9,0.2,0.2, 1e2, 2500., 0.0005);

    // bhorizontal->SetForcingFunction(forcing);    
    //BC
    MatrixDouble val1(3,3);
    val1.setZero();
    VecDouble val2(3);
    val2.setZero();
    //Left bottom
    L2Projection * matbc1 = new L2Projection(3,1,BoundaryConditionType::kDirichlet,val1,val2);
    
    

    cmesh->InsertMaterial(matbc1);

    // val2[0] = 1.e-3;
    val2[1] = -1e1;
    val2[2] = 0.0;
    L2Projection * matbc2 = new L2Projection(4,1,BoundaryConditionType::kNeumann,val1,val2);
    // matbc2->SetForcingFunctionTransient(forcingTransient);
    // matbc2->SetTimeStep(0.5);
    cmesh->InsertMaterial(matbc2);
    

    cmesh->InsertMaterial(bhorizontal);
    // cmesh->InsertMaterial(bvertical);
      
    GmshTools::Read(*cmesh,"../1bar.msh");

    TransientAnalysis an(cmesh,SolverType::ELU,false);
    // LinearAnalysis an(cmesh,SolverType::ELDLt);
    
    std::vector<std::string> ScalarNames, VectorNames;
    VectorNames = {"Displacement"};
    ScalarNames = {"Rotation"};
    an.PrintVariables("result",ScalarNames,VectorNames);
    an.SetMaxIter(10);
    an.Run(500);

    // VTUGenerator::PrintResults(cmesh,"result",ScalarNames,VectorNames);    
        
}           