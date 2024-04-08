
   // Defines the problem dimension
    const int dimension = 2;

{



auto yieldFunction = [](const double &plast, double &sigma_y, double &hardening){
    // sigma_y =.5/sqrt(3.)+500.*plast;
    // sigma_y =50.+500.*plast;
    // hardening = 500.;

    hardening = 121.e3;
    sigma_y = 20000.+hardening*plast;
    // sigma_y = 10.0;//+0.0001*plast;
    
};

    CompMesh* cmesh = new CompMesh(); 

    //Dominio
    // ElasticityPositional2D * matelas = new ElasticityPositional2D(17,20.e6,.25,false);
    Elasticity2D * matelas = new Elasticity2D(17,20.e6,.25,false);
    VonMises *plastmodel = new VonMises(matelas);
    double phi = 40.*M_PI/180.;
    double psi = 00.*M_PI/180.;
    // DruckerPrager *plastmodel = new DruckerPrager(matelas,phi,psi);
    plastmodel->SetUniaxialYieldFunction(yieldFunction);
    // cmesh->InsertMaterial(matelas);
    cmesh->InsertMaterial(plastmodel);

    ElasticTruss *geogrelha = new ElasticTruss(16,2,2.1e7,1);
    cmesh->InsertMaterial(geogrelha);

    //BC
    MatrixDouble val1(2,2);
    val1.setZero();
    VecDouble val2(2);
    VecDouble val3(2);
    val2.setZero();val3.setZero();
    // Fundo 
    // val2[1] = 1.;
    L2Projection * matbc1 = new L2Projection(12,2,0,val1,val2);
    val2.setZero();
    // Laterais
    val2[0] = 1.;
    L2Projection * matbc2 = new L2Projection(13,2,3,val1,val2);
    val2.setZero();
    // Topo
    L2Projection * matbc3 = new L2Projection(14,2,1,val1,val2);
    // Carregamento
    val3[1] = -0.015;
    // val3[1] = -2.e6;
    L2Projection * matbc4 = new L2Projection(15,2,0,val1,val3);

    cmesh->InsertMaterial(matbc1);
    cmesh->InsertMaterial(matbc2);
    cmesh->InsertMaterial(matbc3);
    cmesh->InsertMaterial(matbc4);

    

    // PositionalTruss *geogrid = new PositionalTruss(16, 2, 3000.e12,1.);
    // cmesh->InsertMaterial(geogrid);

    GmshTools::Read(*cmesh,"../geogrelha.msh");
    // GmshTools::Read(*cmesh,"../geogrelha_nosymmetry.msh");
    std::vector<L2Projection *> bcIncrement = {matbc4};

    
    IncrementalAnalysis an(cmesh,SolverType::ELU, 1000, bcIncrement,1.e-6 ,100);
#ifdef RELEASE_BUILD
    an.SType() = SolverType::EUmfpack;
    // an.SType() = SolverType::ECholmod;
#endif



    std::vector<std::string> ScalarNames, VectorNames;
    ScalarNames = {"PlasticStrain"};
    // ScalarNames = {"Pressure","J2","PlasticStrain"};
    // ScalarNames = {"Material"};
    VectorNames = {"Displacement","Stress","RealStress"};
    // VectorNames = {"Displacement"};
    
    an.Run("geogrelha",ScalarNames,VectorNames);

}           