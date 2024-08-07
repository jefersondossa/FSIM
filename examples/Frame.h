
   // Defines the problem dimension
    const int dimension = 1;
    auto forcing = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];
    force[1] = -1;
};

{

    CompMesh* cmesh = new CompMesh(); 

    LinearFrame * bhorizontal = new LinearFrame(5,1.,1.,3.);
    LinearFrame * bvertical = new LinearFrame(4,1.,1.,3.);
    bhorizontal->SetForcingFunction(forcing); 
    //BC
    MatrixDouble val1(3,3);
    val1.setZero();
    VecDouble val2(3);
    val2.setZero();
    //Left bottom
    L2Projection * matbc1 = new L2Projection(3,1,0,val1,val2);

    cmesh->InsertMaterial(matbc1);
    
    cmesh->InsertMaterial(bhorizontal);
    cmesh->InsertMaterial(bvertical);
      
    GmshTools::Read(*cmesh,"../Portico.msh");

    LinearAnalysis an(cmesh,SolverType::ELU);
   
    std::vector<std::string> ScalarNames, VectorNames;
    VectorNames = {"Displacement"};
    ScalarNames = {"Rotation"};
    an.Run();

    VTUGenerator::PrintResults(cmesh,"result",ScalarNames,VectorNames);    
        
}           