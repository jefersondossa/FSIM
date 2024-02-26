
   // Defines the problem dimension
    const int dimension = 2;
{
auto yieldFunction = [](const double &plast, double &sigma_y, double &hardening){
    // sigma_y = 2.e6 - 2.e9*plast;
    sigma_y = 848700;
    // yield = 2.e9;
    hardening = 0.;

};

    CompMesh* cmesh = new CompMesh();

    Elasticity2D * matelas = new Elasticity2D(11,1.e10,.48,false);
    //BC
    MatrixDouble val1(2,2);
    val1.setZero();
    VecDouble val2(2);
    VecDouble val3(2);
    val2.setZero();val3.setZero();
    //Load
    val2[1] = -0.002;
    L2Projection * matbc4 = new L2Projection(6,2,4,val1,val2);
    //Top
    L2Projection * matbc6 = new L2Projection(7,2,3,val1,val2);
    //Left and Right
    val3[0] = 1.;
    L2Projection * matbc3 = new L2Projection(8,2,3,val1,val3);
    L2Projection * matbc5 = new L2Projection(9,2,3,val1,val3);
    //Botom
    val3.setZero();
    val3[1] = 1.;
    L2Projection * matbc1 = new L2Projection(10,2,3,val1,val3);
    

    
    cmesh->InsertMaterial(matbc1);
    cmesh->InsertMaterial(matbc3);
    cmesh->InsertMaterial(matbc4);
    cmesh->InsertMaterial(matbc5);
    cmesh->InsertMaterial(matbc6);
    // MohrCoulomb *plastmodel = new MohrCoulomb(matelas,6.89e4,20.); 
    VonMises *plastmodel = new VonMises(matelas);
    plastmodel->SetUniaxialYieldFunction(yieldFunction);
    // cmesh->InsertMaterial(matelas);
    cmesh->InsertMaterial(plastmodel); 
    std::vector<L2Projection *> bcIncrement = {matbc4};
  
    GmshTools::Read(*cmesh,"../mohrcoulomb.msh");

    IncrementalAnalysis an(cmesh,SolverType::ELDLt, 100, bcIncrement,1.e-5,100);
    // an.SetMaxIter(2);
    std::vector<std::string> ScalarNames, VectorNames;
    ScalarNames = {"PlasticStrain"};
    VectorNames = {"Displacement","Stress","Strain","RealStress"};
    an.Run("plasticitytest",ScalarNames,VectorNames);

    
}           