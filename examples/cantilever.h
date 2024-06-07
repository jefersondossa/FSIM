
   // Defines the problem dimension
    const int dimension = 2;
{
auto yieldFunction = [](const double &plast, double &sigma_y, double &hardening){
    // sigma_y =.5/sqrt(3.)+500.*plast;
    hardening = 0.;
    sigma_y =.24e9+hardening*plast;
    

    // sigma_y = .05;
    // hardening = 0.;
};

    CompMesh* cmesh = new CompMesh(); 

    Elasticity2D * matelas = new Elasticity2D(11,210.e9,.3,true,0.05);
    // ElasticityPositional2D * matelas = new ElasticityPositional2D(10,1.e3,.0,false);

    //BC
    MatrixDouble val1(2,2);
    val1.setZero();
    VecDouble val2(2);
    val2.setZero();
    //Encastre
    val2[0] = 1.;
    L2Projection * matbc1 = new L2Projection(10,2,3,val1,val2);

    val2.setZero();
    //Encastre point
    L2Projection * matbc2 = new L2Projection(8,2,0,val1,val2);

    //Load
    val2[1] = 000.;
    L2Projection * matbc3 = new L2Projection(9,2,1,val1,val2);
    
    cmesh->InsertMaterial(matbc1);
    cmesh->InsertMaterial(matbc2);
    cmesh->InsertMaterial(matbc3);
    
    VonMises *plastmodel = new VonMises(matelas);
    plastmodel->SetUniaxialYieldFunction(yieldFunction);
    // cmesh->InsertMaterial(matelas);
    cmesh->InsertMaterial(plastmodel); 
    std::vector<L2Projection *> bcIncrement = {matbc3};
   
    GmshTools::Read(*cmesh,"../cantilever.msh");

    IncrementalAnalysis an(cmesh,SolverType::ELU, 100, bcIncrement,1.e-7,300);
    VecDouble increment(2);
    // increment.setZero();    
    increment[1] = -1000.;
    // increment[1] = 0.02;
    an.SetIncrement(increment);
#ifdef RELEASE_BUILD
    an.SType() = SolverType::EUmfpack;
    // an.SType() = SolverType::ECholmod;
#endif

    std::vector<std::string> ScalarNames, VectorNames;
    ScalarNames = {"PlasticStrain"};
    VectorNames = {"Displacement","Stress","Strain","RealStress"};
    // VectorNames = {"Displacement","Stress","Strain"};
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