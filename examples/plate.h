
   // Defines the problem dimension
    const int dimension = 2;
{
auto yieldFunction = [](const double &plast, double &sigma_y, double &hardening){
    // sigma_y =.5/sqrt(3.)+500.*plast;
    hardening = 00.;
    sigma_y =16000+hardening*plast;
    

    // sigma_y = .05;
    // hardening = 0.;
};

    CompMesh* cmesh = new CompMesh(); 

    Elasticity2D * matelas = new Elasticity2D(10,1.e7,.24,false);
    //BC
    MatrixDouble val1(2,2);
    val1.setZero();
    VecDouble val2(2);
    val2.setZero();
    //Left
    val2[1] = 0.;
    L2Projection * matbc1 = new L2Projection(9,2,1,val1,val2);

    val2.setZero();
    val2[0] = 0.;
    // //bottom
    // L2Projection * matbc2 = new L2Projection(11,2,1,val1,val2);

    val2.setZero();
    // val2[1] = 0.0005;
    val2[1] = 0.0;
    //Right
    L2Projection * matbc3 = new L2Projection(12,2,1,val1,val2);
    val2.setZero();
    val2[1] = 1.;
    L2Projection * matbc4 = new L2Projection(14,2,3,val1,val2);
    val2.setZero();
    val2[0] = 1.;
    L2Projection * matbc5 = new L2Projection(13,2,3,val1,val2);

    cmesh->InsertMaterial(matbc1);
    // cmesh->InsertMaterial(matbc2);
    cmesh->InsertMaterial(matbc3);
    cmesh->InsertMaterial(matbc4);
    cmesh->InsertMaterial(matbc5);


    VonMises *plastmodel = new VonMises(matelas);
    // Tresca *plastmodel = new Tresca(matelas);
    double angle = 00.*M_PI/180.;
    // DruckerPrager *plastmodel = new DruckerPrager(matelas,angle,angle,false);
    plastmodel->SetUniaxialYieldFunction(yieldFunction);
    // cmesh->InsertMaterial(matelas);
    cmesh->InsertMaterial(plastmodel); 
    std::vector<L2Projection *> bcIncrement = {matbc3};
   
    GmshTools::Read(*cmesh,"../plate.msh");

    IncrementalAnalysis an(cmesh,SolverType::ELU, 50, bcIncrement,1.e-7,30);
    VecDouble increment(2);
    increment.setZero();    
    // increment[1] = 0.00002;
    increment[1] = -5.00;
    an.SetIncrement(increment);
#ifdef RELEASE_BUILD
    an.SType() = SolverType::EUmfpack;
    // an.SType() = SolverType::ECholmod;
#endif

    std::vector<std::string> ScalarNames, VectorNames;
    ScalarNames = {"PlasticStrain"};
    VectorNames = {"Displacement","Stress","Strain"};
    an.Run("plasticitytest",ScalarNames,VectorNames);
    
    // std::vector<std::string> integrate = {"Solution","DerivativeX","DerivativeY"};
    // std::set<int> matid = {6};
    // std::map<std::string,VecDouble> result;
    // cmesh->Integrate(matid,integrate,result);
    // double aux = 206.9/((1.29)*(1.-2.*.29));
    // double valX = aux * ((1.-.29)*result["DerivativeX"][0] + .29*result["DerivativeY"][0])/.45;
    // double valY = aux * (.29*result["DerivativeX"][0] + (1.-.29)*result["DerivativeY"][0])/.45;
    // int a = 0;
    
}           