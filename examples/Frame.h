
   // Defines the problem dimension
    const int dimension = 1;
    auto forcing = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];
    force[1] = 1.;
};

{

    CompMesh* cmesh = new CompMesh(); 

    // LinearFrame * bhorizontal = new LinearFrame(5,1.,1.,3.);
    PositionalFrame2D * bhorizontal = new PositionalFrame2D(5,1.,1.,3.,1.);
    // LinearFrame * bvertical = new LinearFrame(4,1.,1.,3.);
    PositionalFrame2D * bvertical = new PositionalFrame2D(4,1.,1.,3.,1.);
    bvertical->SetForcingFunction(forcing); 
    //BC
    MatrixDouble val1(3,3);
    val1.setZero();
    VecDouble val2(3);
    val2.setZero();
    //Left bottom
    L2Projection * matbc1 = new L2Projection(3,1,BoundaryConditionType::kDirichlet,val1,val2);

    cmesh->InsertMaterial(matbc1);
    
    cmesh->InsertMaterial(bhorizontal);
    cmesh->InsertMaterial(bvertical);
      
    GmshTools::Read(*cmesh,"../Portico.msh");

    NonLinearAnalysis an(cmesh,SolverType::ELDLt,1e-6,1);
   
    std::vector<std::string> ScalarNames, VectorNames;
    VectorNames = {"Displacement"};
    ScalarNames = {"Rotation"};
    an.Run();

    VTUGenerator::PrintResults(cmesh,"result",ScalarNames,VectorNames);    
        
}           