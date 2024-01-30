
   // Defines the problem dimension
    const int dimension = 2;
{
//==========================================================================
//==============================PROBLEM MESHES==============================
//==========================================================================
   
    CompMesh* truss = new CompMesh();
    truss->SetDimension(dimension);
    // PositionalTruss * mattruss = new PositionalTruss(4,dimension,10000.,1.);
    // ElasticTruss * mattruss = new ElasticTruss(4,dimension,1000.,1.);
    // TransientElasticTruss * mattruss = new TransientElasticTruss(4,dimension,1.,1.,0.,1.,0.05);
    TransientPositionalTruss * mattruss = new TransientPositionalTruss(4,dimension,1.,1.,0.,1.0,.05);
    truss->InsertMaterial(mattruss);
    //BC
    MatrixDouble val1(2,2);
    val1.setZero();
    VecDouble val2(2);
    val2.setZero();
    
    L2Projection * matbc1 = new L2Projection(2,2,0,val1,val2);
    val2[0] = 1.;
    L2Projection * matbc2 = new L2Projection(3,2,1,val1,val2);
    
    truss->InsertMaterial(matbc1);
    truss->InsertMaterial(matbc2);
    // std::vector<L2Projection *> bcIncrement = {matbc2};

    // LinearHardening *plastmodel = new LinearHardening(mattruss,111.,4.); 
    // truss->InsertMaterial(plastmodel);
    
    GmshTools::Read(*truss,"../line.msh");
    
    

    
    // LinearAnalysis an(truss,SolverType::ELDLt);
    TransientAnalysis an(truss,SolverType::ELDLt,false);
    // NonLinearAnalysis an(truss,SolverType::ELDLt);
    // NonLinearAnalysis an(arl.MeshVec(),SolverType::EUmfpack,1.e-6,2);
    // NonLinearAnalysis an(coarseModel,SolverType::EUmfpack);
    // an.Run();
    // VecDouble increment(2);
    // increment.setZero();
    // increment[1] = -1.4;
    // an.SetIncrement(increment);

    std::vector<std::string> ScalarNames, VectorNames;
    ScalarNames = {"Stress"};
    VectorNames = {"Displacement"};
    // an.Run();
    an.PrintVariables("truss",ScalarNames,VectorNames);
    an.Run(100);

    // VTUGenerator::PrintResults(truss,"truss",ScalarNames,VectorNames);

}           