
   // Defines the problem dimension
    const int dimension = 2;

{



auto yieldFunction = [](const double &plast, double &sigma_y, double &hardening){

    hardening = 0.;
    sigma_y = 848.7e3+hardening*plast;    
};

    CompMesh* cmesh = new CompMesh(); 

    //Dominio
    ElasticityPositional2D * matelas = new ElasticityPositional2D(12,1.e10,.48,true);
    // Elasticity2D * matelas = new Elasticity2D(12,1.e10,.48,false);
    // VonMises *plastmodel = new VonMises(matelas);
    // double phi = 40.*M_PI/180.;
    // double psi = 00.*M_PI/180.;
    // DruckerPrager *plastmodel = new DruckerPrager(matelas,phi,psi);
    // plastmodel->SetUniaxialYieldFunction(yieldFunction);
    cmesh->InsertMaterial(matelas);
    // cmesh->InsertMaterial(plastmodel);

    //BC
    MatrixDouble val1(2,2);
    val1.setZero();
    VecDouble val2(2);
    VecDouble val3(2);
    val2.setZero();val3.setZero();
    // Fundo 
    // val2[1] = 1.;
    L2Projection * matbc1 = new L2Projection(10,2,0,val1,val2);
    val2.setZero();
    // Laterais
    val2[0] = 1.;
    L2Projection * matbc2 = new L2Projection(9,2,3,val1,val2);
    L2Projection * matbc3 = new L2Projection(11,2,3,val1,val2);
    val2.setZero();
    // Carregamento
    val3[1] = -0.002;
    // val3[1] = -2.e6;
    L2Projection * matbc4 = new L2Projection(8,2,4,val1,val3);

    val2.setZero();
    // Carregamento
    // val3[0] = -1;
    // val3[1] = -2.e6;
    L2Projection * matbc5 = new L2Projection(13,2,1,val1,val3);

    cmesh->InsertMaterial(matbc1);
    cmesh->InsertMaterial(matbc2);
    cmesh->InsertMaterial(matbc3);
    cmesh->InsertMaterial(matbc4);
    cmesh->InsertMaterial(matbc5);

    GmshTools::Read(*cmesh,"../stripfooting.msh");
    std::vector<L2Projection *> bcIncrement = {matbc4};

    

    IncrementalAnalysis an(cmesh,SolverType::ELU, 50, bcIncrement,1.e-6 ,100);
#ifdef RELEASE_BUILD
    an.SType() = SolverType::EUmfpack;
    // an.SType() = SolverType::ECholmod;
#endif

    // VecDouble increment(2);
    // increment.setZero();    
    // // increment[1] = 0.00002;
    // increment[1] = -0.0002;
    // an.SetIncrement(increment);

    std::vector<std::string> ScalarNames, VectorNames;
    // ScalarNames = {"PlasticStrain"};
    // ScalarNames = {"Pressure","J2","PlasticStrain"};
    // ScalarNames = {"Material"};
    // VectorNames = {"Displacement","Stress","RealStress"};
    VectorNames = {"Displacement","Stress"};
    // VectorNames = {"Displacement"};
    
    an.Run("geogrelha",ScalarNames,VectorNames);

}           