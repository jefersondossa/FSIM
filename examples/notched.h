
   // Defines the problem dimension
    const int dimension = 2;
{
auto yieldFunction = [](const REAL &plast, REAL &sigma_y, REAL &hardening){
    // sigma_y = 2.e6 - 2.e9*plast;
    hardening = 0.0;
    sigma_y = 0.45 + hardening*plast;
    // sigma_y = 0.45 + 0.0001*plast;

    //hardening H = d(sigma_y)/d(epsilon_p)
    
    // hardening = 0.0001;
};

    CompMesh* cmesh = new CompMesh();

    Elasticity2D * matelas = new Elasticity2D(11,206.9,.29,false);
    //BC
    MatrixDouble val1(2,2);
    val1.setZero();
    VecDouble val2(2);
    VecDouble val3(2);
    val2.setZero();val3.setZero();
    //Load - Top
    val2[1] = .17;
    L2Projection * matbc4 = new L2Projection(10,2,4,val1,val2);
    //Right and bottom
    L2Projection * matbc6 = new L2Projection(9,2,1,val1,val3);
    L2Projection * matbc3 = new L2Projection(7,2,1,val1,val3);
    //Left
    val3[0] = 1.;
    L2Projection * matbc5 = new L2Projection(8,2,3,val1,val3);
    //Notch
    val3.setZero();
    // val3[1] = 1.;
    L2Projection * matbc1 = new L2Projection(6,2,0,val1,val3);
    

    
    cmesh->InsertMaterial(matbc1);
    cmesh->InsertMaterial(matbc3);
    cmesh->InsertMaterial(matbc4);
    cmesh->InsertMaterial(matbc5);
    cmesh->InsertMaterial(matbc6);
    VonMises *plastmodel = new VonMises(matelas);
    plastmodel->SetUniaxialYieldFunction(yieldFunction);
    // cmesh->InsertMaterial(matelas);
    cmesh->InsertMaterial(plastmodel); 
    std::vector<L2Projection *> bcIncrement = {matbc4};
  
    GmshTools::Read(*cmesh,"../notched.msh");

    IncrementalAnalysis an(cmesh,SolverType::ELDLt, 100, bcIncrement,1.e-6,100  );
#ifdef RELEASE_BUILD
    an.SType() = SolverType::EUmfpack;
    // an.SType() = SolverType::ECholmod;
#endif

    std::vector<std::string> ScalarNames, VectorNames;
    ScalarNames = {"PlasticStrain"};
    VectorNames = {"Displacement","Stress","Strain","RealStress"};
    an.Run("plasticitytest",ScalarNames,VectorNames);
    
    // std::vector<std::string> integrate = {"Solution","DerivativeX","DerivativeY"};
    // std::set<int> matid = {6};
    // std::map<std::string,VecDouble> result;
    // cmesh->Integrate(matid,integrate,result);
    // REAL aux = 206.9/((1.29)*(1.-2.*.29));
    // REAL valX = aux * ((1.-.29)*result["DerivativeX"][0] + .29*result["DerivativeY"][0])/.45;
    // REAL valY = aux * (.29*result["DerivativeX"][0] + (1.-.29)*result["DerivativeY"][0])/.45;
    // int a = 0;
    
}           