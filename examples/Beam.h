
   // Defines the problem dimension
    const int dimension = 1;
    auto forcing = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];
    force[1] = -1;
};

{

    CompMesh* cmesh = new CompMesh(); 

    LinearFrame * matelas = new LinearFrame(4,1.,1.,1.);
    matelas->SetForcingFunction(forcing); 
    //BC
    MatrixDouble val1(3,3);
    val1.setZero();
    VecDouble val2(3);
    val2.setZero();
    //Left
    val2[0] = 1.;
    val2[1] = 1.;
    L2Projection * matbc1 = new L2Projection(2,1,3,val1,val2);

    val2[0] = 0.;
    //Right
    L2Projection * matbc3 = new L2Projection(3,1,3,val1,val2);
    
    cmesh->InsertMaterial(matbc1);
    cmesh->InsertMaterial(matbc3);
 
    cmesh->InsertMaterial(matelas);
      
    GmshTools::Read(*cmesh,"../Beam.msh");

    LinearAnalysis an(cmesh,SolverType::ELU);
   
    std::vector<std::string> ScalarNames, VectorNames;
    VectorNames = {"Displacement"};
    ScalarNames = {"Rotation","BendingMoment","ShearForce"};
    an.Run();

    VTUGenerator::PrintResults(cmesh,"result",ScalarNames,VectorNames);    
        
}           