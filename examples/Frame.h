
   // Defines the problem dimension
    const int dimension = 1;
    auto forcing = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];
    force[1] = 1.;
};

{

    CompMesh* cmesh = new CompMesh(); 

    // LinearFrame * bhorizontal = new LinearFrame(5,1.,1./12.,1.);
    PositionalFrame2D * bhorizontal = new PositionalFrame2D(5,1.,1.,2.);
    // LinearFrame * bvertical = new LinearFrame(4,1.,1.,1.);
    // PositionalFrame2D * bvertical = new PositionalFrame2D(4,1.,1.,1.);
    // bvertical->SetForcingFunction(forcing); 
    //BC
    MatrixDouble val1(3,3);
    val1.setZero();
    VecDouble val2(3);
    val2.setZero();
    //Left bottom
    L2Projection * matbc1 = new L2Projection(3,1,BoundaryConditionType::kDirichlet,val1,val2);

    cmesh->InsertMaterial(matbc1);

    val2[1] = 1.;
    L2Projection * matbc2 = new L2Projection(4,1,BoundaryConditionType::kNeumann,val1,val2);

    cmesh->InsertMaterial(matbc2);
    
    cmesh->InsertMaterial(bhorizontal);
    // cmesh->InsertMaterial(bvertical);
      
    GmshTools::Read(*cmesh,"../1bar.msh");

    NonLinearAnalysis an(cmesh,SolverType::ELU,1e-6,2);
    // LinearAnalysis an(cmesh,SolverType::ELDLt);
   
    std::vector<std::string> ScalarNames, VectorNames;
    VectorNames = {"Displacement"};
    ScalarNames = {"Rotation"};
    an.Run();

    VTUGenerator::PrintResults(cmesh,"result",ScalarNames,VectorNames);    
        
}           