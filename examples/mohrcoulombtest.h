
   // Defines the problem dimension
    const int dimension = 2;
{
auto yieldFunction = [](const double &plast, double &sigma_y){
    sigma_y = 2.e6 - 2.e9*plast;
};

    CompMesh* cmesh = new CompMesh();

    Elasticity2D * matelas = new Elasticity2D(11,20.e9,.2);
    //BC
    MatrixDouble val1(2,2);
    val1.setZero();
    VecDouble val2(2);
    VecDouble val3(2);
    val2.setZero();val3.setZero();
    val2[0] = 0.001;
    L2Projection * matbc4 = new L2Projection(9,1,0,val1,val2);
    
    L2Projection * matbc6 = new L2Projection(12,1,0,val1,val3);

    val3[0] = 1.;
    L2Projection * matbc3 = new L2Projection(8,1,3,val1,val3);
    // val3[0] = 1.;
    L2Projection * matbc5 = new L2Projection(10,1,1,val1,val3);
    L2Projection * matbc1 = new L2Projection(6,1,1,val1,val3);
    

    
    cmesh->InsertMaterial(matbc1);
    cmesh->InsertMaterial(matbc3);
    cmesh->InsertMaterial(matbc4);
    cmesh->InsertMaterial(matbc5);
    cmesh->InsertMaterial(matbc6);
    // MohrCoulomb *plastmodel = new MohrCoulomb(matelas,6.89e4,20.); 
    // LinearHardening *plastmodel = new LinearHardening(matelas,6.89e4,1.e6); 
    VonMises *plastmodel = new VonMises(matelas);
    plastmodel->SetUniaxialYieldFunction(yieldFunction);
    // cmesh->InsertMaterial(matelas);
    cmesh->InsertMaterial(plastmodel); 
    std::vector<L2Projection *> bcIncrement = {matbc4};
  
    GmshTools::Read(*cmesh,"../mohrcoulomb.msh");

    IncrementalAnalysis an(cmesh,SolverType::ELDLt, 20, bcIncrement,1.e-4);
    
    std::vector<std::string> ScalarNames, VectorNames;
    ScalarNames = {};
    VectorNames = {"Displacement","Stress","Strain"};
    an.Run("plasticitytest",ScalarNames,VectorNames);

    
}           