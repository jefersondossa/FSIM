
   // Defines the problem dimension
    const int dimension = 2;

{

auto yieldFunction = [](const double &plast, double &sigma_y, double &hardening){
    // sigma_y =.5/sqrt(3.)+500.*plast;
    sigma_y =50.+500.*plast;
    hardening = 500.;

    // sigma_y = .05;
    // hardening = 0.;
};

    CompMesh* cmesh = new CompMesh(); 

    //Dominio
    // ElasticityPositional2D * matelas = new ElasticityPositional2D(17,1000.,.0);
    Elasticity2D * matelas = new Elasticity2D(17,1000.,.0,false);
    // VonMises *plastmodel = new VonMises(matelas);
    DruckerPrager *plastmodel = new DruckerPrager(matelas,0.,0.);
    plastmodel->SetUniaxialYieldFunction(yieldFunction);
    // cmesh->InsertMaterial(matelas);
    // LinearHardening *plastmodel = new LinearHardening(matelas,111.,4.); 
    cmesh->InsertMaterial(plastmodel);

    //BC
    MatrixDouble val1(2,2);
    val1.setZero();
    VecDouble val2(2);
    VecDouble val3(2);
    val2.setZero();val3.setZero();
    // Fundo 
    L2Projection * matbc1 = new L2Projection(12,2,0,val1,val2);
    // Laterais
    L2Projection * matbc2 = new L2Projection(13,2,0,val1,val2);
    // Topo
    L2Projection * matbc3 = new L2Projection(14,2,1,val1,val2);
    // Carregamento
    val3[1] = -300.;
    L2Projection * matbc4 = new L2Projection(15,2,1,val1,val3);

    cmesh->InsertMaterial(matbc1);
    cmesh->InsertMaterial(matbc2);
    cmesh->InsertMaterial(matbc3);
    cmesh->InsertMaterial(matbc4);

    

    PositionalTruss *geogrid = new PositionalTruss(16, 2, 3000.,1.);
    cmesh->InsertMaterial(geogrid);

    GmshTools::Read(*cmesh,"../geogrelha.msh");
    std::vector<L2Projection *> bcIncrement = {matbc4};

    
    IncrementalAnalysis an(cmesh,SolverType::ELU, 30, bcIncrement);
    
    std::vector<std::string> ScalarNames, VectorNames;
    ScalarNames = {"PlasticStrain"};
    VectorNames = {"Displacement","Stress","Strain"};
    
    an.Run("geogrelha",ScalarNames,VectorNames);

}           